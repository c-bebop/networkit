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
    AlgebraicAPD(const Graph& graph) : At(Matrix::adjacencyMatrix(graph, MinPlusSemiring::zero()).transpose())
    {
        
    }

    void run()
    {

    }

private:
    Graph g;
    const Matrix At;
};

} /* namespace NetworKit */