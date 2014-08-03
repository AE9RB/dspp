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

namespace testing {
    
    /// @class Test benchtest.hpp
    /// @ingroup benchtest
    /// @brief Abstract class to be inherited by test fixtures.
    /// @details See \ref benchtest "benchtest module documentation"
    /// for more usage information.
    /// ~~~
    /// class MyCase : public testing::Test {
    /// protected:
    ///     void MyTest() {
    ///     }
    /// }
    /// TEST_F(MyCase, MyTest)
    /// ~~~
    class Test {
        friend class Runner;
        virtual Info* TestInfo() = 0;
        virtual void TestBody() {};
        // Helper for common typos.
        struct Setup_should_be_spelled_SetUp {};
        virtual Setup_should_be_spelled_SetUp* Setup() {
            return nullptr;
        }
        virtual Setup_should_be_spelled_SetUp* SetupTestCase() {
            return nullptr;
        }
        // Benchmark info
        std::chrono::high_resolution_clock::time_point start_time;
        ::std::vector<double> results;
        unsigned long count = 0;
        unsigned long size = 0;
        bool keep_running = true;
    protected:
        Test() {}
        /// @ingroup benchtest
        /// @brief Called once before any tests in this case are run.
        static void SetUpTestCase() {}
        /// @ingroup benchtest
        /// @brief Called once when all tests in this case are finished.
        static void TearDownTestCase() {}
        /// @ingroup benchtest
        /// @brief Called at the start of every test using this fixture.
        virtual void SetUp() {}
        /// @ingroup benchtest
        /// @brief Called at the end of every test using this fixture.
        virtual void TearDown() {}
        /// @ingroup benchtest
        /// @brief Controls execution of a benchmarked section.
        /// @param max Limit the maximum number of executions.
        /// @retval bool True indicates the benchmarked section should execute again.
        bool Benchmark(unsigned int max = 100) {
            auto end_time = ::std::chrono::high_resolution_clock::now();
            double this_time = ::std::chrono::duration<double, std::micro>(end_time - start_time).count();
            if (!keep_running) return false;
            if (max < 5) max = 5;
            if (!size) {
                size = max / 5;
                if (size < 5) size = 5;
            }
            if (count) {
                if (results.size() < size) {
                    auto it=results.begin();
                    for (; it!=results.end(); ++it) {
                        if (*it > this_time) {
                            results.insert(it, this_time);
                            break;
                        }
                    }
                    if (it==results.end()) {
                        results.push_back(this_time);
                    }
                }
                else {
                    keep_running = false;
                    for (auto it=results.begin(); it!=results.end(); ++it) {
                        if (*it > this_time) {
                            results.pop_back();
                            results.insert(it, this_time);
                            keep_running = true;
                            break;
                        }
                    }
                }
            }
            ++count;
            if (count > max) keep_running = false;
            if (!keep_running) {
                double mean = 0;
                unsigned long trim = results.size() / 5;
                for (auto i = trim; i < results.size() - trim; ++i) {
                    mean += results[i];
                }
                mean /= results.size() - trim*2;
                reporter()->Bench(count-1, mean);
            }
            start_time = ::std::chrono::high_resolution_clock::now();
            return keep_running;
        }
    public:
        /// True if a fatal failure (ASSERT) happened in this test.
        bool HasFatalFailure() {
            return TestInfo()->HasFatalFailure();
        }
        /// True if a non-fatal failure (EXPECT) happened in this test.
        bool HasNonfatalFailure() {
            return TestInfo()->HasNonfatalFailure();
        }
        /// True if either kind of failure was encountered.
        bool HasFailure() {
            return TestInfo()->HasFailure();
        }
    };

}
