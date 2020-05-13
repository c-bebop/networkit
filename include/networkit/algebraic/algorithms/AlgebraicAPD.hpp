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

    }

private:
    const Matrix matrix;
};

} /* namespace NetworKit */