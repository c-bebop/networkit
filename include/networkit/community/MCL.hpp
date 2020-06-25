#pragma once

#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>

#include <vector>
#include <stdexcept>
#include <cmath>

namespace NetworKit {

template<class Matrix>
Vector columnSum(Matrix const& m)
{
    Vector w(m.numberOfColumns());

    for (size_t i = 0; i < m.numberOfRows(); ++i)
    {
        m.forElementsInRow(i, [&](index column, double value) 
        {
            w[column] += value;
        });
    }

    return w;
}

template<class Matrix>
void normalize(Matrix& m)
{
    Vector w = columnSum<Matrix>(m);
    w.forElements([](double& x) { x = 1. / x; });

    for (size_t i = 0; i < m.numberOfRows(); ++i)
    {
        m.forNonZeroElementsInRow(i, [&](index j, double value)
        {
            m.setValue(i, j, value * (0.01 / w[j]));
        });
    }
}

class MCL : public CommunityDetectionAlgorithm {

public:
    MCL(const Graph& G, double delta = 0.01, size_t expansion = 3, double inflation = 1.5)
    : CommunityDetectionAlgorithm(G)
    , m_delta(delta)
    , m_expansion(expansion)
    , m_inflation(inflation)
    {
        if (m_expansion < 2)
        {
            throw std::runtime_error("Expansion factor can not be scmaller 2.");
        }

        if (m_inflation <= 1.)
        {
            throw std::runtime_error("Inflaction factor must be greater 1.");
        }
    }
    
    template<class Matrix>
    bool equals(Matrix const& a, Matrix const& b, double delta)
    {
        bool recurse = false;

        for (size_t i = 0; i < a.numberOfRows() && !recurse; ++i)
        {
            a.forNonZeroElementsInRow(i, [&](index j, double value)
            {
                double const epsilon = std::abs(a(i, j) - b(i, j));
                if (epsilon > delta)
                {
                    recurse = true;
                }
            });
        }

        return recurse;
    }

    template <typename Matrix>
    void inflation(Matrix& m, double r)
    {
        for (size_t i = 0; i < m.numberOfRows(); ++i)
        {
            m.forElementsInRow(i, [&](index j, double value) {
                m.setValue(i, j, value * r);
            });
        }
    }

    template<typename Matrix>
    Matrix mcl(Matrix const& M, size_t e, double r)
    {
        bool recurse = false;

        Matrix C_i = M;
        Matrix C_f;
        do
        {
            // Actual algorithm
            // expansion
            C_f = C_i * C_i;
            for (size_t i = 1; i < e; ++i)
            {
                C_f = C_f * C_f;
            }

            inflation(C_f, m_inflation);


            // Normalise Columns
            Vector w = columnSum<Matrix>(C_f);
            w.forElements([](double& x) { x = 1. / x; });

            for (size_t i = 0; i < C_f.numberOfRows(); ++i)
            {
                C_f.forNonZeroElementsInRow(i, [&](index j, double value)
                {
                    C_f.setValue(i, j, C_f(i, j) * w[j]);
                });
            }

            recurse = equals(C_f, C_i, m_delta);

            if (recurse)
            {
                C_i = C_f;
            }

        } while (recurse);

        return std::move(C_f);
    }

    void run() override
    {
        // Creating Markov Matrix
        CSRMatrix A = CSRMatrix::adjacencyMatrix(*G);

        // Transpose because A(i, j) is the weight from j to i to generate the Markov Matrix
        // Maybe not needed? Is that only a detail?!
        // A.transpose();
        
        Vector v = columnSum<CSRMatrix>(A);

        CSRMatrix D = CSRMatrix::diagonalMatrix(v);
        // invserse D
        for (size_t i = 0; i < D.numberOfRows(); ++i)
        {
            D.forElementsInRow(i, [&](index j, double value)
            {
                D.setValue(i, j, 1. / D(i, j)); 
            });
        }

        CSRMatrix M_sparse = A * D;
        DenseMatrix M(M_sparse.numberOfRows(), M_sparse.numberOfColumns());
        for (size_t i = 0; i < M_sparse.numberOfRows(); ++i)
        {
            M_sparse.forElementsInRow(i, [&](index j, double value) {
                M.setValue(i, j, M_sparse(i, j));
            });
        }

        DenseMatrix C_f = mcl(M, m_expansion, m_delta);

        index z = G->upperNodeIdBound();
        if (result.numberOfElements() != z) {
            result = Partition(z);
            result.allToSingletons();
        }

        for (size_t i = 0; i < C_f.numberOfRows(); ++i)
        {
            if (C_f(i, i) > m_delta)
            {
                std::cout << "C_f(i, i) > m_delta" << std::endl;
                C_f.forElementsInRow(i, [&](size_t j, double x)
                {
                    if (x > m_delta && i != j)
                    {
                        result.moveToSubset(j, i);
                    }
                });
            }
        }

        hasRun = true;
    }

    std::string toString() const override { return std::string("MCL"); }

    /**
    * Get number of iterations in last run.
    *
    * @return The number of iterations.
    */
    count numberOfIterations() { return 0; };

    /**
    * Get list of running times for each iteration.
    *
    * @return The list of running times in milliseconds
    */
    std::vector<count> getTiming() { return std::vector<count>{1, 2, 3, 4}; }


private:
    double m_delta{0.01};
    size_t m_expansion{2};
    double m_inflation{1.01};
};

}
