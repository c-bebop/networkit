#include <gtest/gtest.h>

#include <networkit/algebraic/DynamicMatrix.hpp>

#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/algebraic/algorithms/AlgebraicAPD.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>

namespace NetworKit {

class AlgebraicAPDGTest : public testing::Test {
public:
    AlgebraicAPDGTest() = default;
    virtual ~AlgebraicAPDGTest() = default;

protected:
    std::vector<double> classicBF(const Graph& graph, node s) const;
};

TEST_F(AlgebraicAPDGTest, first_test) {
    KONECTGraphReader reader;
    Graph graph = reader.read("/Users/vivid/HU/GALA/networkit_fork_rep/input/out.moreno_lesmis_lesmis");
    std::cout << "graph has " << graph.numberOfNodes() << " nodes and " << graph.numberOfEdges() 
              << " edges and directed? " << graph.isDirected() << std::endl;

    node const n = graph.numberOfNodes();

    std::vector<double> entries(n * n, 0.);
    
    for (node i = 0; i < n; ++i)
    {
        for (node j = 0; j < n; ++j)
        {
            entries[i * n + j] = graph.hasEdge(i, j) ? 1. : 0.;
        }
    }

    DenseMatrix dense_matrix(n, n, entries);

    AlgebraicAPD<DenseMatrix> apd(dense_matrix);
}

} /* namespace NetworKit */
