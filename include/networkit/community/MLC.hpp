#pragma once

#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>

#include <vector>
#include <stdexcept>
#include <cmath>

namespace NetworKit {

class MCL : public CommunityDetectionAlgorithm {

public:
    MCL(const Graph& G, double delta = 0.1, size_t expansion = 2, double inflation = 1.01)
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
        

    void run() override
    {

        // Creating Markov Matrix
        CSRMatrix A = CSRMatrix::adjacencyMatrix(*G);

        // Transpose becouse A(i, j) is the weight from j to i to generate the Markov Matrix
        // Maybe not needed? Is that only a detail?!
        A.transpose();
        
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

        // Actual algorithm
        // expansion
        DenseMatrix C_f = M * M;
        for (size_t i = 1; i < m_expansion; ++i)
        {
            C_f = C_f * C_f;
        }

        // inflation
        for (size_t i = 0; i < C_f.numberOfRows(); ++i)
        {
            C_f.forElementsInRow(i, [&](index j, double value) {
                C_f.setValue(i, j, C_f(i, j) * m_inflation);
            });
        }

        Vector w = columnSum<DenseMatrix>(C_f);

        // Normalise Columns
        w /= 1.;
        for (size_t i = 0; i < C_f.numberOfRows(); ++i)
        {
            C_f.forNonZeroElementsInRow(i, [&](index j, double value)
            {
                C_f.setValue(i, j, C_f(i, j) * w[j]);
            });
        }

        bool recurse = false;
        for (size_t i = 0; i < C_f.numberOfRows() && !recurse; ++i)
        {
            C_f.forNonZeroElementsInRow(i, [&](index j, double value)
            {
                if (!recurse)
                {
                    double const epsilon = std::abs(C_f(i, j) - M(i, j));
                    recurse = epsilon > m_delta;
                }
            });
        }
    }

    std::string toString() const override { return std::string("MLC"); }

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
    double m_delta{0.1};
    size_t m_expansion{2};
    double m_inflation{1.01};
};

}
