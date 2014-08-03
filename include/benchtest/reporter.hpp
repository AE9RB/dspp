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

    //TODO enable argv check for xcode in Error()

    /// @class Reporter benchtest.hpp
    /// @brief Abstract class for building custom reporters.
    /// @details A custom reporter might be desired for XML or JSON output
    /// to be parsed by another tool. You may also build a custom reporter
    /// by inheriting DefaultReporter if you just want a few adjustments.
    class Reporter {
    public:
        /// An instance of Info for the current test is put here
        /// automatically by the framework.
        class Info* test_info;
        /// Start is called before any tests are run.
        virtual void Start(size_t cases, size_t total_qty) = 0;
        /// End is called after all tests are run.
        virtual void End(long ms) = 0;
        /// StartCase is called before any tests in a case are run.
        virtual void StartCase(size_t case_qty) = 0;
        /// EndCase is called after all tests in a case are finished.
        virtual void EndCase(long ms) = 0;
        /// Run is called once before each test.
        virtual void Run() = 0;
        /// Pass is called at the end of a test that didn't report failures.
        virtual void Pass(long ms) = 0;
        /// Fail is called at the end of a test that failed.
        virtual void Fail(long ms) = 0;
        /// Bench is called during test to report a benchmark time.
        virtual void Bench(long iterations, double us) = 0;
        /// Print is called with miscellaneous annotations for a human-readable report.
        virtual void Print(::std::string message) = 0;
        /// Trace is called to print each line in a scoped trace.
        virtual void Trace(::std::string message, const char* file, long line) = 0;
        /// Error is called for each failure during a test.
        virtual void Error(::std::string message, const char* file, long line) = 0;
    };

    /// @class DefaultReporter benchtest.hpp
    /// @brief Reporter class that outputs a human readable format.
    class DefaultReporter : public Reporter {
    public:
        /// Creates an instance which can optionally report somewhere other than STDOUT.
        DefaultReporter(::std::ostream& ostream = ::std::cout) {
            this->ostream = &ostream;
        }
    protected:
        /// The ostream this instance was initialized to.
        ::std::ostream* ostream;
        /// Remember the number of cases passed to Start().
        size_t cases;
        /// Remember the total_qty of tests passed to Start().
        size_t total_qty;
        /// Remember the number of tests for the current test case.
        size_t case_qty;
        /// Remember the name of every failed test for the final report.
        ::std::vector<::std::string> failures;

        /// Utility for adding "s", e.g., "case" becomes "cases" when not qty 1.
        virtual ::std::string Pluralize(size_t qty, const char* label = nullptr) {
            auto str = ::std::to_string(qty);
            str += " ";
            if (label != nullptr) str += label;
            else str += "test";
            if (qty != 1) {
                if (str.back() < 'a') str += "S";
                else str += "s";
            }
            return str;
        }

        virtual void Start(size_t cases, size_t total_qty) {
            this->cases = cases;
            this->total_qty = total_qty;
            failures.clear();
            *ostream << "[==========] Running " << Pluralize(total_qty) << " from ";
            *ostream << Pluralize(cases, "test case") << "." << ::std::endl;
        }

        virtual void End(long ms) {
            *ostream << "[==========] " << Pluralize(total_qty) << " from ";
            *ostream << Pluralize(cases, "test case") << " ran. ";
            *ostream << "(" << ms << " ms total)" << ::std::endl;
            *ostream << "[  PASSED  ] " << Pluralize(total_qty - failures.size()) << "." << ::std::endl;
            if (!failures.empty()) {
                *ostream << "[  FAILED  ] " << Pluralize(failures.size()) << ", listed below:" << ::std::endl;
                for (auto f : failures) {
                    *ostream << "[  FAILED  ] " << f << ::std::endl;
                }
                *ostream << ::std::endl << " " << Pluralize(failures.size(), "FAILED TEST") << ::std::endl;
            }
        }

        virtual void StartCase(size_t case_qty) {
            this->case_qty = case_qty;
            *ostream << "[----------] " << Pluralize(case_qty) << " from " << test_info->test_case_name() << ::std::endl;

        }

        virtual void EndCase(long ms) {
            *ostream << "[----------] " << Pluralize(case_qty) << " from " << test_info->test_case_name();
            *ostream << " (" << ms << " ms total)" << ::std::endl << ::std::endl;
        }

        virtual void Run() {
            *ostream << "[ RUN      ] " << test_info->name() << ::std::endl;
        }

        virtual void Pass(long ms) {
            *ostream << "[       OK ] " << test_info->name() << " (" << ms << " ms)" << ::std::endl;
        }

        virtual void Fail(long ms) {
            failures.push_back(test_info->name());
            *ostream << "[  FAILED  ] " << test_info->name() << " (" << ms << " ms)" << ::std::endl;
        }
        
        virtual void Bench(long iterations, double us) {
            *ostream << "[   TIME   ] " << iterations << " iterations, " << us << " us" << ::std::endl;
        }

        virtual void Print(::std::string message) {
            *ostream << message << ::std::endl;
        }

        virtual void Trace(::std::string message, const char* file, long line) {
            *ostream << file << ":" << line << ": " << message << ::std::endl;
        }

        virtual void Error(::std::string message, const char* file, long line) {
            *ostream << file << ":" << line << ": ";
            if (false) {
                // For Apple Xcode, you can run script in build phase:
                // "${TARGET_BUILD_DIR}/${EXECUTABLE_NAME}" --xcode
                // Benchmarks will be skipped and failures will
                // be annotated like compiler errors.
                for (auto ch : message) {
                    *ostream << ch;
                    if (ch == '\n') *ostream << file << ":" << line << ": ";
                }
                *ostream << ::std::endl;
            } else {
                *ostream << "Failure" << ::std::endl;
                *ostream << message << ::std::endl;
            }
        }

    };

}
