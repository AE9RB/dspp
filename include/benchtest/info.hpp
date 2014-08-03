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
    
    /// @class Info benchtest.hpp
    /// @brief Information about current test for the Reporter.
    /// @details An instance of this will automatically be attached to
    /// the Reporter. It is used to obtain strings for the report output.
    class Info {
        friend class Runner;
        const char* case_name;
        const char* test_name;
        const char* type_name;
        virtual class Test* CreateFixture() = 0;
    protected:
        /// @cond BENCHTEST_IMPL
        Info(const char* case_name, const char* test_name, const char* type_name) :
        case_name(case_name),
        test_name(test_name),
        type_name(type_name) {
            if (::std::string(type_name).size() == 0) {
                this->type_name = nullptr;
            }
        }
    public:
        virtual void SetUpTestCase() = 0;
        virtual void TearDownTestCase() = 0;
        int fatal_failure_count = 0;
        int nonfatal_failure_count = 0;
        bool HasFatalFailure() {
            return fatal_failure_count;
        }
        bool HasNonfatalFailure() {
            return nonfatal_failure_count;
        }
        bool HasFailure() {
            return HasFatalFailure() || HasNonfatalFailure();
        }
        /// @endcond

        /// Name of the test case including template type.
        ::std::string test_case_name() const {
            auto n = ::std::string(case_name);
            if (type_name) {
                n += "<";
                n += type_name;
                n += ">";
            }
            return n;
        }

        /// Full name of the test, i.e., case_name + type_name + test_name.
        ::std::string name() const {
            auto n = test_case_name();
            n += ".";
            n += test_name;
            return n;
        }
        
    };
    
    
}