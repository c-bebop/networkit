#pragma once

#include <cassert>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/algebraic/GraphBLAS.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Implementation of the Bellman-Ford algorithm using the GraphBLAS interface.
 */
template<class Matrix>
class AlgebraicAPD : public Algorithm {
public:
    AlgebraicAPD(const Matrix& m)
    : matrix(m)
    {
        
    }

    void run()
    {
        distance = apd(matrix);
    }

    Matrix the_new_a(Matrix a, Matrix z)
    {
        std::vector<double> entries(a.numberOfRows() * a.numberOfColumns(), 0.);
        for (count i = 0; i < a.numberOfRows(); ++i)
        {
            for (count j = 0; j < a.numberOfColumns(); ++j)
            {
                if (i != j && (a(i, j) == 1 || z(i, j) > 0))
                {
                    entries[i * a.numberOfRows() + j] = 1;
                }
            }
        }
        
        return Matrix(a.numberOfRows(), a.numberOfColumns(), entries);
    }

    Matrix apd(Matrix a)
    {
        Matrix z = a * a;
        Matrix a_ = the_new_a(a, z);

        return a_;
    }

    Matrix distance;
private:
    const Matrix matrix;
};

} /* namespace NetworKit */