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

    class Reporter {
    public:
        class Info* test_info;
        virtual void Start(size_t cases, size_t total_qty) = 0;
        virtual void End(long ms) = 0;
        virtual void StartCase(size_t case_qty) = 0;
        virtual void EndCase(long ms) = 0;
        virtual void Run() = 0;
        virtual void Pass(long ms) = 0;
        virtual void Fail(long ms) = 0;
        virtual void Print(::std::string message) = 0;
        virtual void Trace(::std::string message, const char* file, long line) = 0;
        virtual void Error(::std::string message, const char* file, long line) = 0;
    };


    class DefaultReporter : public Reporter {
    public:
        DefaultReporter(::std::ostream& ostream = ::std::cout) {
            this->ostream = &ostream;
        }
    private:
        ::std::ostream* ostream;
        size_t cases;
        size_t total_qty;
        size_t case_qty;
        ::std::vector<::std::string> failures;

        ::std::string Pluralize(size_t qty, const char* label = nullptr) {
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

        void Start(size_t cases, size_t total_qty) {
            this->cases = cases;
            this->total_qty = total_qty;
            failures.clear();
            *ostream << "[==========] Running " << Pluralize(total_qty) << " from ";
            *ostream << Pluralize(cases, "test case") << "." << ::std::endl;
        }

        void End(long ms) {
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

        void StartCase(size_t case_qty) {
            this->case_qty = case_qty;
            *ostream << "[----------] " << Pluralize(case_qty) << " from " << test_info->test_case_name() << ::std::endl;

        }

        void EndCase(long ms) {
            *ostream << "[----------] " << Pluralize(case_qty) << " from " << test_info->test_case_name();
            *ostream << " (" << ms << " ms total)" << ::std::endl << ::std::endl;
        }

        void Run() {
            *ostream << "[ RUN      ] " << test_info->name() << ::std::endl;
        }

        void Pass(long ms) {
            *ostream << "[       OK ] " << test_info->name() << " (" << ms << " ms)" << ::std::endl;
        }

        void Fail(long ms) {
            failures.push_back(test_info->name());
            *ostream << "[  FAILED  ] " << test_info->name() << " (" << ms << " ms)" << ::std::endl;
        }

        void Print(::std::string message) {
            *ostream << message << ::std::endl;
        }

        void Trace(::std::string message, const char* file, long line) {
            *ostream << file << ":" << line << ": " << message << ::std::endl;
        }

        void Error(::std::string message, const char* file, long line) {
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
