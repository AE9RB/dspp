// benchtest - A benchmarking and unit testing framework.
// Copyright (C) 2014 David Turnbull
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "benchtest.hpp"

// Self-test for benchtest.
// There is no Makefile. Simply run this self-test like so:
// g++ -Wall -o selftest -std=c++11 -O3 selftest.cpp && ./selftest
// Since failures must be tested, it is normal see them in the report.
// Message to STDERR and exit code are used to determine sucess.

#define SELFTEST_ABORT \
{std::cerr << std::endl << "ERROR: self-test failed.\n"; abort();}

#define SELFTEST_BASE1 \
do { \
bool selftest_died_ = true; \
auto selftest_lambda_ = [&]() {

#define SELFTEST_BASE2(fatal, qty) \
;selftest_died_ = false;}; \
auto selftest_fails1_ = ::testing::reporter()->test_info->fatal_failure_count + \
::testing::reporter()->test_info->nonfatal_failure_count; \
selftest_lambda_(); \
auto selftest_fails2_ = ::testing::reporter()->test_info->fatal_failure_count + \
::testing::reporter()->test_info->nonfatal_failure_count; \
if(selftest_died_ != fatal || selftest_fails2_ - selftest_fails1_ != qty) \
SELFTEST_ABORT \
} while (false);

#define SELFTEST_FATAL(expr) \
SELFTEST_BASE1 expr SELFTEST_BASE2(true, 1)

#define SELFTEST_NONFATAL(expr) \
SELFTEST_BASE1 expr SELFTEST_BASE2(false, 1)

#define SELFTEST_FATAL_QTY2(expr) \
SELFTEST_BASE1 expr SELFTEST_BASE2(true, 2)

#define SELFTEST_NONFATAL_QTY2(expr) \
SELFTEST_BASE1 expr SELFTEST_BASE2(false, 2)

#define SELFTEST_SUCCESS(expr) \
SELFTEST_BASE1 expr SELFTEST_BASE2(false, 0)


int main() {
    testing::reporter(new testing::DefaultReporter);
    auto result = testing::Runner::RunAll();
    if (result == EXIT_SUCCESS) SELFTEST_ABORT;
    std::cerr << std::endl << "Self-test was successful." << std::endl;
    return 0; // Clean exit code if no problems.
}

// Program flow and control tests

TEST(Control, TraceAndExplicit) {
    {
        SCOPED_TRACE();
        SELFTEST_FATAL( FAIL() << "Game on!"; );
        SCOPED_TRACE() << "Whoa" << ", cowboy!";
        SELFTEST_NONFATAL( ADD_FAILURE(); );
    }
    SELFTEST_NONFATAL( ADD_FAILURE_AT("file.bogus", -1); );
}

void no_fatal(bool assert_true) {
    EXPECT_TRUE(true);
    ASSERT_TRUE(assert_true);
}

#include <cmath>
TEST(Control, NoFatalFailure) {
    SELFTEST_SUCCESS(
        ASSERT_NO_FATAL_FAILURE(no_fatal(true));
    )
    SELFTEST_SUCCESS(
        EXPECT_NO_FATAL_FAILURE(no_fatal(true));
    )
    SELFTEST_FATAL_QTY2(
        ASSERT_NO_FATAL_FAILURE(no_fatal(false));
    )
    SELFTEST_NONFATAL_QTY2(
        EXPECT_NO_FATAL_FAILURE(no_fatal(false));
    )
}

// Exceptions

TEST(Exceptions, Throw) {
    SELFTEST_SUCCESS(
        EXPECT_THROW(throw std::runtime_error("foo"), std::runtime_error);
    )
    SELFTEST_SUCCESS(
        ASSERT_THROW(throw std::runtime_error("foo"), std::runtime_error);
    )
    SELFTEST_NONFATAL(
        EXPECT_THROW({}, std::runtime_error);
    )
    SELFTEST_FATAL(
        ASSERT_THROW({}, std::runtime_error);
    )
    SELFTEST_NONFATAL(
        EXPECT_THROW(throw std::exception(), std::runtime_error);
    )
    SELFTEST_FATAL(
        ASSERT_THROW(throw std::exception(), std::runtime_error);
    )
}

TEST(Exceptions, AnyThrow) {
    SELFTEST_SUCCESS(
        EXPECT_ANY_THROW(throw std::exception());
    )
    SELFTEST_SUCCESS(
        ASSERT_ANY_THROW(throw std::runtime_error("foo"));
    )
    SELFTEST_NONFATAL(
        EXPECT_ANY_THROW({});
    )
    SELFTEST_FATAL(
        ASSERT_ANY_THROW({});
    )
}

TEST(Exceptions, NoThrow) {
    SELFTEST_SUCCESS(
        EXPECT_NO_THROW({});
    )
    SELFTEST_SUCCESS(
        ASSERT_NO_THROW({});
    )
    SELFTEST_NONFATAL(
        EXPECT_NO_THROW(throw std::exception());
    )
    SELFTEST_FATAL(
        ASSERT_NO_THROW(throw std::runtime_error("foo"));
    )
}

// Fixture tests

class Fixture1 : public testing::Test {
protected:
    static int test_case_counter;
    static void SetUpTestCase() {
        ADD_FAILURE() << "Non-fatal failure in SetUpTestCase";
        ++test_case_counter;
    }
    int test_counter = 0;
    void SetUp() {
        ADD_FAILURE() << "Non-fatal failure in SetUp";
        ++test_counter;
    }
};
int Fixture1::test_case_counter = 0;

TEST_F(Fixture1, NonFatalSetup) {
    if (1 != test_case_counter) SELFTEST_ABORT;
    if (1 != test_counter) SELFTEST_ABORT;
    //    ADD_FAILURE() << "TestBody called";
}

class Fixture2 : public testing::Test {
protected:
    static void TearDownTestCase() {
        FAIL();
    }
    void TearDown() {
        FAIL();
    }
    void SetUp() {
        FAIL();
    }
};

TEST_F(Fixture2, FatalSetup) {
    SELFTEST_ABORT;
}

class Fixture3 : public testing::Test {
protected:
    static void SetUpTestCase() {
        FAIL();
    }
};

TEST_F(Fixture3, FatalSetUpTestCase) {
    SELFTEST_ABORT;
}


template<typename T>
class TypedFixture : public testing::Test {
    static T dat;
    T val = 0;
protected:
    static void SetUpTestCase() {
        dat = 7;
    }
    void SetUp() {
        val = 9;
    }
    void Test() {
        SELFTEST_SUCCESS( EXPECT_EQ(7, dat) )
        SELFTEST_SUCCESS( EXPECT_EQ(9, val) )
    }
};
template<typename T>
T TypedFixture<T>::dat = 0;

TEST_T(TypedFixture, int, Test)
TEST_T(TypedFixture, float, Test);

// Predecate tests

bool pred1(bool v1) {
    return v1;
}
bool pred2(bool v1, int v2) {
    return v1;
}
bool pred3(bool v1, int v2, float v3) {
    return v1;
}

TEST(Predecate, Bool) {
    SELFTEST_SUCCESS(
        EXPECT_PRED1(pred1, true);
        EXPECT_PRED2(pred2, true, 2);
        EXPECT_PRED3(pred3, true, 2, 3.0);
        ASSERT_PRED1(pred1, true);
        ASSERT_PRED2(pred2, true, 2);
        ASSERT_PRED3(pred3, true, 2, 3.0);
    )
    SELFTEST_NONFATAL( EXPECT_PRED1(pred1, false); );
    SELFTEST_NONFATAL( EXPECT_PRED2(pred2, false, 2); );
    SELFTEST_NONFATAL( EXPECT_PRED3(pred3, false, 2, 3.0); );
    SELFTEST_FATAL( ASSERT_PRED1(pred1, false); );
    SELFTEST_FATAL( ASSERT_PRED2(pred2, false, 2); );
    SELFTEST_FATAL( ASSERT_PRED3(pred3, !true, 2+1, 3.0001); );
}

// Boolean tests

TEST(TrueFalse, Passing) {
    SELFTEST_SUCCESS(
        EXPECT_TRUE(1==1);
        EXPECT_FALSE(1==0);
        ASSERT_TRUE(1==1);
        ASSERT_FALSE(1==0);
    )
}

TEST(TrueFalse, Failing) {
    SELFTEST_NONFATAL( EXPECT_TRUE(1==0); );
    SELFTEST_NONFATAL( EXPECT_FALSE(1==1); );
    SELFTEST_FATAL( ASSERT_TRUE(1==0); );
    SELFTEST_FATAL( ASSERT_FALSE(1==1); );
}

::testing::AssertionResult IsEven(int n) {
    if ((n % 2) == 0)
        return ::testing::AssertionSuccess() << n << " is even";
    else
        return ::testing::AssertionFailure() << n << " is odd";
}

TEST(TrueFalse, AssertionResult) {
    SELFTEST_SUCCESS(
        EXPECT_FALSE(IsEven(1));
        EXPECT_TRUE(IsEven(2));
    )
    SELFTEST_NONFATAL( EXPECT_TRUE(IsEven(3)); );
    SELFTEST_NONFATAL( EXPECT_FALSE(IsEven(4)); );
    SELFTEST_FATAL( ASSERT_TRUE(IsEven(5)); );
    SELFTEST_FATAL( ASSERT_FALSE(IsEven(6)); );

}

// Comparison tests

TEST(Comparisons, Equality) {
    SELFTEST_SUCCESS(
        EXPECT_EQ(std::string("X"), std::string("X"));
        ASSERT_EQ(4.5, 4.1 + 0.4);
        EXPECT_NE(5, 3+3);
        ASSERT_NE(9.9, 3+3);
    )
    SELFTEST_NONFATAL( EXPECT_EQ(std::string("X"), std::string("Y")); )
    SELFTEST_FATAL( ASSERT_EQ(9.9, 3+3); )
    SELFTEST_NONFATAL( EXPECT_NE(5.0, 4 + 1); )
    SELFTEST_FATAL( ASSERT_NE(4.5, 4.1 + 0.4); )
}

TEST(Comparisons, GreaterLess) {
    SELFTEST_SUCCESS(
        EXPECT_GT(5.1, 4 + 1);
        ASSERT_GT(4.6, 4.1 + 0.4);
        EXPECT_GE(5.0, 4 + 1);
        ASSERT_GE(4.5, 4.1 + 0.4);
        EXPECT_LT(5, 3+3);
        ASSERT_LT(9.9, 7+3);
        EXPECT_LE(5, 3+3);
        ASSERT_LE(6, 3+3);
    )
    SELFTEST_NONFATAL( EXPECT_GT(5, 3+3); )
    SELFTEST_FATAL( ASSERT_GT(3.9, 3+3); )
    SELFTEST_NONFATAL( EXPECT_GE(4.9, 4 + 1); )
    SELFTEST_FATAL( ASSERT_GE(4.4, 4.1 + 0.4); )
    SELFTEST_NONFATAL( EXPECT_LT(7, 3+3); )
    SELFTEST_FATAL( ASSERT_LT(9.9, 3+3); )
    SELFTEST_NONFATAL( EXPECT_LE(5.1, 4 + 1); )
    SELFTEST_FATAL( ASSERT_LE(4.6, 4.1 + 0.4); )
}

TEST(Comparisons, Floats) {
    SELFTEST_NONFATAL( EXPECT_EQ(float(0.5), double(0.50000000001)); )
    SELFTEST_SUCCESS( EXPECT_EQ(0.5f, 0.5000001f); )
    SELFTEST_NONFATAL( EXPECT_EQ(0.5f, 0.500001f); )
    SELFTEST_FATAL( ASSERT_EQ(0.5f, 0.500001f); )
    SELFTEST_SUCCESS( EXPECT_EQ(-0.5, -0.5000000000000001); )
    SELFTEST_NONFATAL( EXPECT_EQ(-0.5, -0.500000000000001); )
    SELFTEST_FATAL( ASSERT_EQ(-0.5, -0.500000000000001); )
}

TEST(Comparisons, Near) {
    SELFTEST_SUCCESS( EXPECT_NEAR(500, 499, 1); )
    SELFTEST_NONFATAL( EXPECT_NEAR(500+1, 500-1, 1); )
    SELFTEST_NONFATAL( EXPECT_NEAR(500, 500.0-0.6, 0.01); )
    SELFTEST_FATAL( ASSERT_NEAR(500.0+0.6, 500.0-0.6, 1); )
}
