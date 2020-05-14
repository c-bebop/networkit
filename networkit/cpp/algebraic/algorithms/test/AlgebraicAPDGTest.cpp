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

TEST_F(AlgebraicAPDGTest, homework) 
{
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

TEST_F(AlgebraicAPDGTest, small_example) 
{
    Graph graph(5);
    graph.addEdge(0, 1);
    graph.addEdge(0, 3);
    graph.addEdge(1, 2);
    graph.addEdge(1, 3);
    graph.addEdge(4, 2);

    AlgebraicAPD apd(graph);
    apd.run();

    DenseMatrix distance = apd.getDistance();

    EXPECT_EQ(0, distance(0, 0));
    EXPECT_EQ(1, distance(0, 1));
    EXPECT_EQ(2, distance(0, 2));
    EXPECT_EQ(1, distance(0, 3));
    EXPECT_EQ(3, distance(0, 4));

    EXPECT_EQ(1, distance(1, 0));
    EXPECT_EQ(0, distance(1, 1));
    EXPECT_EQ(1, distance(1, 2));
    EXPECT_EQ(1, distance(1, 3));
    EXPECT_EQ(2, distance(1, 4));
    
    EXPECT_EQ(2, distance(2, 0));
    EXPECT_EQ(1, distance(2, 1));
    EXPECT_EQ(0, distance(2, 2));
    EXPECT_EQ(2, distance(2, 3));
    EXPECT_EQ(1, distance(2, 4));
    
    EXPECT_EQ(1, distance(3, 0));
    EXPECT_EQ(1, distance(3, 1));
    EXPECT_EQ(2, distance(3, 2));
    EXPECT_EQ(0, distance(3, 3));
    EXPECT_EQ(3, distance(3, 4));
    
    EXPECT_EQ(3, distance(4, 0));
    EXPECT_EQ(2, distance(4, 1));
    EXPECT_EQ(1, distance(4, 2));
    EXPECT_EQ(3, distance(4, 3));
    EXPECT_EQ(0, distance(4, 4));
}

} /* namespace NetworKit */
