#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/algebraic/MatrixTools.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>

#include <chrono>
#include <networkit/algebraic/SPA.hpp>
#include <memory>

namespace NetworKit {

inline void printmatrix(CSRMatrix const& m)
{
    for (int i = 0; i < m.numberOfRows(); ++i)
    {
        for (int j = 0; j < m.numberOfColumns(); ++j)
        {
            std::cout << m(i, j);
        }

        std::cout << std::endl;
    }
}

class CSRMatrixSpGEMMGTest : public testing::Test {
public:
    CSRMatrixSpGEMMGTest() = default;
    virtual ~CSRMatrixSpGEMMGTest() = default;
};


TEST_F(CSRMatrixSpGEMMGTest, SPA_reset)
{
    SPA spa(10);
    
    spa.accumulate(1, 1.);

    ASSERT_EQ(1u, spa.nnz());

    spa.reset();

    ASSERT_EQ(0u, spa.nnz());
}

TEST_F(CSRMatrixSpGEMMGTest, SPA_output)
{
    std::vector<Triplet> triplets = {{0,0,1}, {0,1,2}, {0,2,3}, {1,0,2}, {1,1,2}, {2,0,3}, {2,2,3}, {2,3,-1}, {3,2,-1}, {3,3,4}};

    //
    //				 1  2  3  0
    // 				 2  2  0  0
    // mat1 = mat2 = 3  0  3 -1
    //				 0  0 -1  4
    //
    CSRMatrix expected(4, triplets);
    ASSERT_EQ(4u, expected.numberOfRows());
    ASSERT_EQ(4u, expected.numberOfColumns());

    SPA spa(expected.numberOfColumns());

    std::vector<index> rowIdx(expected.numberOfRows() + 1, 0);
    std::vector<index> columnIdx;
    std::vector<double> nonZeros;

    //				 1  2  3  0
    spa.accumulate(1., 0);
    spa.accumulate(2., 1);
    spa.accumulate(3., 2);

    size_t nnz = spa.output(nonZeros, columnIdx, 0);
    EXPECT_EQ(3, nnz);
    size_t current_nnz = nnz;
    rowIdx[1] = current_nnz;
    spa.reset();

    // 				 2  2  0  0
    spa.accumulate(2., 0);
    spa.accumulate(2., 1);

    nnz = spa.output(nonZeros, columnIdx, nonZeros.size());
    EXPECT_EQ(2, nnz);
    current_nnz += nnz;
    rowIdx[2] = current_nnz;
    spa.reset();

    //               3  0  3 -1
    spa.accumulate(3., 0);
    spa.accumulate(3., 2);
    spa.accumulate(-1., 3);

    nnz = spa.output(nonZeros, columnIdx, nonZeros.size());
    EXPECT_EQ(3, nnz);
    current_nnz += nnz;
    rowIdx[3] = current_nnz;
    spa.reset();
    
    //				 0  0 -1  4
    spa.accumulate(-1., 2);
    spa.accumulate(4., 3);

    nnz = spa.output(nonZeros, columnIdx, nonZeros.size());
    EXPECT_EQ(2, nnz);
    current_nnz += nnz;
    rowIdx[4] = current_nnz;
    spa.reset();
    
    CSRMatrix result(4, 4, rowIdx, columnIdx, nonZeros);

    EXPECT_TRUE(expected == result);
}

TEST_F(CSRMatrixSpGEMMGTest, SpGEMM_SPA_minimal_example)
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


TEST_F(CSRMatrixSpGEMMGTest, SpGEMM_SPA_airfoil)
{
    METISGraphReader reader;
    Graph G = reader.read("input/airfoil1.graph");

    CSRMatrix A = CSRMatrix::adjacencyMatrix(G);

    auto spa_timer_begin = std::chrono::high_resolution_clock::now();
    CSRMatrix AA = A.spgemm_spa(A);
    auto spa_timer_end = std::chrono::high_resolution_clock::now();

    std::cout << "SPGEMM SPA time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(spa_timer_end - spa_timer_begin).count() << " ns" << std::endl;

    auto mm_timer_begin = std::chrono::high_resolution_clock::now();
    CSRMatrix AA_mm = A * A;
    auto mm_timer_end = std::chrono::high_resolution_clock::now();

    std::cout << "MM time:         " << std::chrono::duration_cast<std::chrono::nanoseconds>(mm_timer_end - mm_timer_begin).count() << " ns" << std::endl;

    EXPECT_TRUE(AA_mm == AA);
}


} /* namespace NetworKit */
