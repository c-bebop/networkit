#pragma once

#include <cassert>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>

namespace NetworKit {

/**
 * @ingroup algebraic
 * Implementation of the Bellman-Ford algorithm using the GraphBLAS interface.
 */
class AlgebraicAPD : public Algorithm {
public:
    AlgebraicAPD(const Graph& graph)
    {
        node const n = graph.numberOfNodes();

        std::vector<double> entries(n * n, 0.);
        
        for (node i = 0; i < n; ++i)
        {
            for (node j = 0; j < n; ++j)
            {
                entries[i * n + j] = graph.hasEdge(i, j) ? 1. : 0.;
            }
        }

        matrix = DenseMatrix(n, n, entries);
    }

    void run()
    {
        distance = apd(matrix);
    }

    DenseMatrix the_new_a(DenseMatrix const& a, DenseMatrix const& z)
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
        
        return DenseMatrix(a.numberOfRows(), a.numberOfColumns(), entries);
    }

    bool all_nodes_reached(DenseMatrix const& a)
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

    DenseMatrix make_distance(DenseMatrix const& d_, DenseMatrix const& s, DenseMatrix const& z)
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

        return DenseMatrix(d_.numberOfRows(), d_.numberOfColumns(), entries);
    }

    DenseMatrix apd(DenseMatrix const& a)
    {
        DenseMatrix const z = a * a;
        DenseMatrix const a_ = the_new_a(a, z);

        if (all_nodes_reached(a_))
        {
            return (a_ * 2) - a;
        }

        DenseMatrix const d_ = apd(a_);
        DenseMatrix const s = a * d_;

        return make_distance(d_, s, z);
    }

    DenseMatrix distance;
private:
    DenseMatrix matrix;
};

} /* namespace NetworKit */