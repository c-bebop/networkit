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

    Matrix the_new_a(Matrix const& a, Matrix const& z)
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

    bool all_nodes_reached(Matrix const& a)
    {
        bool reached = true;
        for (count i = 0; i < a.numberOfRows() && reached; ++i)
        {
            for (count j = 0; j < a.numberOfColumns() && reached; ++j)
            {
                if (i != j)
                {
                    reached = a(i, j) == 1;
                }
            }
        }
        return reached;
    }

    Matrix make_distance(Matrix const& d_, Matrix const& s, Matrix const& z)
    {
        std::vector<double> entries(d_.numberOfRows() * d_.numberOfColumns(), 0.);
        for (count i = 0; i < d_.numberOfRows(); ++i)
        {
            for (count j = 0; j < d_.numberOfColumns(); ++j)
            {
                if (i != j)
                {
                    entries[i * d_.numberOfRows() + j] = [&]() {
                        if (s(i, j) >= (d_(i, j) * z(i, i)))
                        {
                            return 2 * d_(i, j);
                        }

                        if (s(i, j) < (d_(i, j) * z(i, i)))
                        {
                            return 2 * d_(i, j) - 1.;
                        }

                        return 0.;
                    }();
                }
            }
        }

        return Matrix(d_.numberOfRows(), d_.numberOfColumns(), entries);
    }

    Matrix apd(Matrix const& a)
    {
        Matrix const z = a * a;
        Matrix const a_ = the_new_a(a, z);

        if (all_nodes_reached(a_))
        {
            return (a_ * 2) - a;
        }

        Matrix const d_ = apd(a_);
        Matrix const s = a * d_;

        return make_distance(d_, s, z);
    }

    Matrix distance;
private:
    const Matrix matrix;
};

} /* namespace NetworKit */