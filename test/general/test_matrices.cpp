#include "gtest/gtest.h"
#include "cantera/numerics/BandMatrix.h"
#include "cantera/numerics/DenseMatrix.h"

using namespace Cantera;

class BandMatrixTest : public testing::Test
{
public:
    BandMatrixTest()
        : x(vector_fp{1,2,3,4,5,6})
        , b1(vector_fp{-8, -8, 6, 40, 149, 81})
        , b2(vector_fp{1, 8, 30, 72, 140, 65})
        , v1(vector_fp{3, 13, 28, 55, 100, 92})
        , v2(vector_fp{0, 5, 16, 39, 115, 116})
    {
        A1.resize(6, 1, 2); // one lower, two upper
        A2.resize(6, 2, 1); // two lower, one upper
        // main diagonal
        for (int i = 0; i < 6; i++) {
            A1(i, i) = i + 1;
            A2(i, i) = i + 1;
        }

        // first subdiagonal and superdiagonal
        for (int i = 0; i < 5; i++) {
            A1(i+1, i) = 2 * i + 1;
            A1(i, i+1) = i * i;
            A2(i+1, i) = 2 * i + 1;
            A2(i, i+1) = i * i;
        }

        // second subdiagonal and superdiagonal
        for (int i = 0; i < 4; i++) {
            A1(i, i+2) = - i - 3;
            A2(i+2, i) = - i - 1;
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

class DenseMatrixTest : public testing::Test
{
public:
    DenseMatrixTest()
        : x4(vector_fp{1,2,3,4})
        , x3(vector_fp{3,2,1})
        , b1(vector_fp{14, 32, 66, -34})
        , b2(vector_fp{-6, 4, 26, 2})
        , b3(vector_fp{14, 32, 66})
    {
        A1.resize(4, 4); // square
        A2.resize(4, 3); // more rows
        A3.resize(3, 4); // more columns

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                A1(i,j) = A2(i,j) = A3(i,j) = i*i + 2*j*j + i*j - 3;
            }
            A1(3,i) = A2(3,i) = pow(-1, i);
            A1(i,3) = A3(i,3) = pow(2, i);
        }
        A1(3,3) = -9;
    }

    double special_sum(DenseMatrix M) {
        double sum = 0;
        for (size_t i = 0; i < M.nRows(); i++) {
            for (size_t j = 0; j < M.nColumns(); j++) {
                sum += M(i,j) * (i+2*j+1);
            }
        }
        return sum;
    }

    DenseMatrix A1, A2, A3;
    vector_fp x4, x3;
    vector_fp b1, b2, b3;
};

TEST_F(DenseMatrixTest, matrix_times_vector)
{
    vector_fp c(4, 0.0);
    A1.mult(x4.data(), c.data());
    for (size_t i = 0; i < 4; i++) {
        EXPECT_DOUBLE_EQ(b1[i], c[i]);
    }
    A2.mult(x3.data(), c.data());
    for (size_t i = 0; i < 4; i++) {
        EXPECT_DOUBLE_EQ(b2[i], c[i]);
    }
    A3.mult(x4.data(), c.data());
    for (size_t i = 0; i < 3; i++) {
        EXPECT_DOUBLE_EQ(b3[i], c[i]);
    }
}

TEST_F(DenseMatrixTest, matrix_times_matrix)
{
    DenseMatrix c(4, 3);
    A1.mult(A2, c);
    EXPECT_DOUBLE_EQ(3033, special_sum(c));

    c.resize(3, 4);
    A3.mult(A1, c);
    EXPECT_DOUBLE_EQ(3386, special_sum(c));

    c.resize(3, 3);
    A3.mult(A2, c);
    EXPECT_DOUBLE_EQ(2989, special_sum(c));

    c.resize(4, 4);
    A2.mult(A3, c);
    EXPECT_DOUBLE_EQ(4014, special_sum(c));
}

TEST_F(DenseMatrixTest, solve_single_rhs)
{
    vector_fp c(b1);
    solve(A1, c.data());
    for (size_t i = 0; i < 4; i++) {
        EXPECT_NEAR(x4[i], c[i], 1e-12);
    }
}

TEST_F(DenseMatrixTest, solve_multi_rhs)
{
    DenseMatrix B(A1.nColumns(), 5);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            B(i,j) = b1[i] * (j+1);
        }
    }
    solve(A1, B);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 5; j++) {
            EXPECT_NEAR(x4[i] * (j+1), B(i,j), 1e-12);
        }
    }
}

TEST_F(DenseMatrixTest, increment)
{
    vector_fp c(b1.size(), 3.0);
    increment(A1, x4.data(), c.data());
    for (size_t i = 0; i < 4; i++) {
        EXPECT_DOUBLE_EQ(3.0 + b1[i], c[i]);
    }

    c.assign(b2.size(), 3.0);
    increment(A2, x3.data(), c.data());
    for (size_t i = 0; i < 4; i++) {
        EXPECT_DOUBLE_EQ(3.0 + b2[i], c[i]);
    }

    c.assign(b3.size(), 3.0);
    increment(A3, x4.data(), c.data());
    for (size_t i = 0; i < 3; i++) {
        EXPECT_DOUBLE_EQ(3.0 + b3[i], c[i]);
    }
}

TEST_F(DenseMatrixTest, invert_full)
{
    DenseMatrix B(A1);
    DenseMatrix C(A1.nRows(), A1.nColumns());
    invert(B);
    A1.mult(B, C);
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            if (i == j) {
                EXPECT_NEAR(1.0, C(i,j), 1e-14);
            } else {
                EXPECT_NEAR(0.0, C(i,j), 1e-14);
            }
        }
    }
}

TEST_F(DenseMatrixTest, invert_partial)
{
    DenseMatrix B(A1);
    DenseMatrix Aref(A1);
    size_t N = 3;
    invert(B, N);

    DenseMatrix As(3, 3);
    DenseMatrix Bs(3, 3);
    DenseMatrix C(3, 3);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            As(i,j) = A1(i,j);
            Bs(i,j) = B(i,j);
        }
    }

    As.mult(Bs, C);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            if (i == j) {
                EXPECT_NEAR(1.0, C(i,j), 1e-14);
            } else {
                EXPECT_NEAR(0.0, C(i,j), 1e-14);
            }
        }
    }
    for (size_t i = 0; i < 4; i++) {
        EXPECT_DOUBLE_EQ(Aref(3,i), A1(3,i));
        EXPECT_DOUBLE_EQ(Aref(i,3), A1(i,3));
    }
}
