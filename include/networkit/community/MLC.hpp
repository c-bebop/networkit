#pragma once

#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>

#include <vector>

namespace NetworKit {

class MCL : public CommunityDetectionAlgorithm {

public:
    MCL(const Graph& G, size_t expansion = 2, double inflation = 0.1, count delta = none)
    : CommunityDetectionAlgorithm(G)
    , m_delta(delta)
    , m_expansion(expansion)
    , m_inflation(inflation)
    {

    }

    void run() override
    {
        // Creating Markov Matrix
        CSRMatrix A = CSRMatrix::adjacencyMatrix(*G);

        // Transpose becouse A(i, j) is the weight from j to i to generate the Markov Matrix
        A.transpose();
        
        Vector v(A.numberOfColumns());

        for (size_t i = 0; i < A.numberOfRows(); ++i)
        {
            A.forElementsInRow(i, [&](index column, double value) 
            {
                v[column] += value;
            });
        }

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

        // expansion
        for (size_t i = 0; i < m_expansion; ++i)
        {
            M = M * M;
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
    count m_delta{0};
    size_t m_expansion{2};
    double m_inflation{0.1};
};

}
