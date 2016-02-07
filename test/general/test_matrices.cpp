#include "gtest/gtest.h"
#include "cantera/numerics/BandMatrix.h"

using namespace Cantera;

class BandMatrixTest : public testing::Test
{
public:
    BandMatrixTest()
        : x{1,2,3,4,5,6}
        , b1{-8, -8, 6, 40, 149, 81}
        , b2{1, 8, 30, 72, 140, 65}
        , v1{3, 13, 28, 55, 100, 92}
        , v2{0, 5, 16, 39, 115, 116}
    {
        A1.resize(6, 1, 2); // one lower, two upper
        A2.resize(6, 2, 1); // two lower, one upper
        // main diagonal
        for (size_t i = 0; i < 6; i++) {
            A1(i, i) = i + 1;
            A2(i, i) = i + 1;
        }

        // first subdiagonal and superdiagonal
        for (size_t i = 0; i < 5; i++) {
            A1(i+1, i) = 2 * i + 1;
            A1(i, i+1) = i * i;
            A2(i+1, i) = 2 * i + 1;
            A2(i, i+1) = i * i;
        }

        // second subdiagonal and superdiagonal
        for (size_t i = 0; i < 4; i++) {
            A1(i, i+2) = - static_cast<int>(i + 3);
            A2(i+2, i) = - static_cast<int>(i + 1);
        }
    }

    BandMatrix A1, A2;
    vector_fp x;
    vector_fp b1, b2;
    vector_fp v1, v2;
};

TEST_F(BandMatrixTest, matrix_times_vector)
{
    vector_fp c(6, 0.0);
    A1.mult(x.data(), c.data());
    for (size_t i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(b1[i], c[i]);
    }
    A2.mult(x.data(), c.data());
    for (size_t i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(b2[i], c[i]);
    }
}

TEST_F(BandMatrixTest, vector_times_matrix)
{
    vector_fp c(6, 0.0);
    A1.leftMult(x.data(), c.data());
    for (size_t i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(v1[i], c[i]);
    }
    A2.leftMult(x.data(), c.data());
    for (size_t i = 0; i < 6; i++) {
        EXPECT_DOUBLE_EQ(v2[i], c[i]);
    }
}

TEST_F(BandMatrixTest, solve_linear_system)
{
    vector_fp c(6, 0.0);
    A1.solve(b1.data(), c.data());
    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR(x[i], c[i], 1e-10);
    }
    A2.solve(b2.data(), c.data());
    for (size_t i = 0; i < 6; i++) {
        EXPECT_NEAR(x[i], c[i], 1e-10);
    }
}

TEST_F(BandMatrixTest, oneNorm) {

    EXPECT_DOUBLE_EQ(28, A1.oneNorm());
    EXPECT_DOUBLE_EQ(23, A2.oneNorm());
}

TEST_F(BandMatrixTest, checkRowsColumns) {
    double s;
    size_t i = A1.checkRows(s);
    EXPECT_EQ((size_t) 0, i);
    EXPECT_DOUBLE_EQ(3, s);
    i = A2.checkRows(s);
    EXPECT_EQ((size_t) 0, i);
    EXPECT_DOUBLE_EQ(1, s);

    i = A1.checkColumns(s);
    EXPECT_EQ((size_t) 0, i);
    EXPECT_DOUBLE_EQ(1, s);
    i = A2.checkColumns(s);
    EXPECT_EQ((size_t) 0, i);
    EXPECT_DOUBLE_EQ(1, s);
}
