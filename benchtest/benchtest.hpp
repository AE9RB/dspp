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

#ifndef BENCHTEST_HPP
#define BENCHTEST_HPP

#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <complex>
#include <cmath>

namespace testing {

    template<bool=true>
    class Reporter* reporter(class Reporter *newr = nullptr) {
        static class Reporter* r;
        if (newr != nullptr) r = newr;
        return r;
    }

}

#include "info.hpp"
#include "printer.hpp"
#include "reporter.hpp"
#include "asserter.hpp"
#include "test.hpp"
#include "runner.hpp"


#define BENCHTEST_CLASS_NAME_(test_case_name, test_name, type_name) \
test_case_name##_##test_name##_TestBench_##type_name

#define BENCHTEST_(test_case_name, test_name, parent_class, type_name) \
friend class ::testing::Factory<BENCHTEST_CLASS_NAME_(test_case_name, test_name, type_name)>; \
static ::testing::Info* const test_info; \
::testing::Info* TestInfo() { return this->test_info;} \
void TestBody();\
};\
::testing::Info* const BENCHTEST_CLASS_NAME_(test_case_name, test_name, type_name)::test_info = \
new ::testing::Factory<BENCHTEST_CLASS_NAME_(test_case_name, test_name, type_name)>(#test_case_name, #test_name, #type_name); \
void BENCHTEST_CLASS_NAME_(test_case_name, test_name, type_name)::TestBody()

#define BENCHTEST_F_(test_case_name, test_name, parent_class) \
class BENCHTEST_CLASS_NAME_(test_case_name, test_name, ) : public parent_class { \
BENCHTEST_(test_case_name, test_name, parent_class, )

#define BENCHTEST_T_(test_case_name, test_name, parent_class, type_name) \
class BENCHTEST_CLASS_NAME_(test_case_name, test_name, type_name) : public parent_class<type_name> { \
BENCHTEST_(test_case_name, test_name, parent_class, type_name)

#define TEST(test_case_name, test_name) \
BENCHTEST_F_(test_case_name, test_name, ::testing::Test)

#define TEST_F(test_fixture, test_name) \
BENCHTEST_F_(test_fixture, test_name, test_fixture)

#define TEST_T(test_fixture, type_name, method_name) \
BENCHTEST_T_(test_fixture, method_name, test_fixture, type_name) {method_name();}


#endif
