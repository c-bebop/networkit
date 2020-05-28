#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/algebraic/MatrixTools.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>

#include <networkit/algebraic/SPA.hpp>
#include <memory>

namespace NetworKit {

class CSRMatrixSpGEMMGTest : public testing::Test {
public:
    CSRMatrixSpGEMMGTest() = default;
    virtual ~CSRMatrixSpGEMMGTest() = default;
};


TEST_F(CSRMatrixSpGEMMGTest, SPA_reset)
{
    SPAGala spa(10);
    
    spa.accumulate(1, 1.);

    ASSERT_EQ(1u, spa.nnz());

    spa.reset();

    ASSERT_EQ(0u, spa.nnz());
}
TEST_F(CSRMatrixSpGEMMGTest, SPA_output)
{
    SPAGala spa(10);
    
    spa.accumulate(1, 1.);
    spa.accumulate(3, 3.);
    spa.accumulate(5, 5.);

    EXPECT_EQ(3u, spa.nnz());

    std::vector<std::pair<size_t, double>> row = spa.output();

    EXPECT_EQ(row[0].first, 1); 
    EXPECT_EQ(row[0].second, 1.); 

    EXPECT_EQ(row[1].first, 3); 
    EXPECT_EQ(row[1].second, 3.); 

    EXPECT_EQ(row[2].first, 5); 
    EXPECT_EQ(row[2].second, 5.); 
}

TEST_F(CSRMatrixSpGEMMGTest, test)
{
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    //
    //				 1  2  3  0
    // 				 2  2  0  0
    // mat1 = mat2 = 3  0  3 -1
    //				 0  0 -1  4
    //
    CSRMatrix mat1(4, triplets);
    ASSERT_EQ(4u, mat1.numberOfRows());
    ASSERT_EQ(4u, mat1.numberOfColumns());

    CSRMatrix mat2(4, triplets);
    ASSERT_EQ(4u, mat2.numberOfRows());
    ASSERT_EQ(4u, mat2.numberOfColumns());

    //
    //			14  6  12  -3
    //			 6  8   6   0
    // result = 12  6  19  -7
    //			-3  0  -7  17
    //
    CSRMatrix result = mat1.spgemm_spa(mat2);
    ASSERT_EQ(mat1.numberOfRows(), result.numberOfRows());
    ASSERT_EQ(mat1.numberOfColumns(), result.numberOfColumns());

    EXPECT_EQ(14, result(0,0));
    EXPECT_EQ(6, result(0,1));
    EXPECT_EQ(12, result(0,2));
    EXPECT_EQ(-3, result(0,3));
    EXPECT_EQ(6, result(1,0));
    EXPECT_EQ(8, result(1,1));
    EXPECT_EQ(6, result(1,2));
    EXPECT_EQ(0, result(1,3));
    EXPECT_EQ(12, result(2,0));
    EXPECT_EQ(6, result(2,1));
    EXPECT_EQ(19, result(2,2));
    EXPECT_EQ(-7, result(2,3));
    EXPECT_EQ(-3, result(3,0));
    EXPECT_EQ(0, result(3,1));
    EXPECT_EQ(-7, result(3,2));
    EXPECT_EQ(17, result(3,3));
}


} /* namespace NetworKit */
