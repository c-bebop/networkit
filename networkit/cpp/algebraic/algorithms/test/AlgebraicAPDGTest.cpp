#include <gtest/gtest.h>

#include <networkit/distance/Dijkstra.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>

#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/algebraic/algorithms/AlgebraicAPD.hpp>


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
    Graph graph = reader.read("/Users/vivid/HU/GALA/moreno_lesmis/out.moreno_lesmis_lesmis");

    AlgebraicAPD<CSRMatrix> bf(graph);
}

} /* namespace NetworKit */
