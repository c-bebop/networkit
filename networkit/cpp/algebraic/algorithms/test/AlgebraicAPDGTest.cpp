#include <gtest/gtest.h>

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
    Graph graph = reader.read("input/out.moreno_lesmis_lesmis");
    std::cout << "graph has " << graph.numberOfNodes() << " nodes and " << graph.numberOfEdges() 
              << " edges and directed? " << graph.isDirected() << std::endl;

    AlgebraicAPD apd(graph);
    apd.run();

    DenseMatrix distance = apd.getDistance();

    for (count i = 0; i < distance.numberOfRows(); ++i)
    {
        std::cout << "Distance[" << i + 1 << "]: ";
        for (count j = 0; j < distance.numberOfColumns(); ++j)
        {
            std::cout << distance(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

} /* namespace NetworKit */
