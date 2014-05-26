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

    //TODO validate all cases really use same fixture
    //TODO add argv options
    
    class Runner {
        
        template <typename T>
        static int RunTime(::std::chrono::time_point<T> start_time) {
            auto end_time = ::std::chrono::high_resolution_clock::now();
            return ::std::chrono::duration<double, ::std::milli>(end_time - start_time).count();
        }
        
        static ::std::map<::std::string, ::std::list<class Info*>>& testers() {
            static ::std::map<::std::string, ::std::list<class Info*>> testers;
            return testers;
        }

    public:
        
        static void AddTest(Info* tester) {
            testers()[tester->test_case_name()].push_back(tester);
        }

#if defined(__GNUC__) && !defined(COMPILER_ICC)
        __attribute__ ((warn_unused_result))
#endif
        static int RunAll()
        {
            auto total_start = ::std::chrono::high_resolution_clock::now();
            bool has_failures = false;
            size_t numtests = 0;
            for (auto& x : testers()) {
                numtests += x.second.size();
            }
            reporter()->Start(testers().size(), numtests);
            for (auto& tester : testers()) {
                auto case_start = ::std::chrono::high_resolution_clock::now();
                reporter()->test_info = tester.second.front();
                reporter()->StartCase(tester.second.size());
                reporter()->test_info->SetUpTestCase();
                has_failures = has_failures || reporter()->test_info->HasFailure();
                bool set_up_case_fatal_failure = reporter()->test_info->HasFatalFailure();
                if (set_up_case_fatal_failure) {
                    ::std::string s("Fatal failure in ");
                    s += tester.first;
                    s += "::SetUpTestCase(); not running tests.";
                    reporter()->Print(s);
                }
                for (auto& test_info : tester.second) {
                    auto test_start = ::std::chrono::high_resolution_clock::now();
                    reporter()->test_info = test_info;
                    if (set_up_case_fatal_failure) {
                        reporter()->Fail(RunTime(test_start));
                    } else {
                        reporter()->Run();
                        auto test = test_info->CreateFixture();
                        test->SetUp();
                        if (!reporter()->test_info->HasFatalFailure()) test->TestBody();
                        test->TearDown();
                        if (reporter()->test_info->HasFailure()) {
                            has_failures = true;
                            reporter()->Fail(RunTime(test_start));
                        } else {
                            reporter()->Pass(RunTime(test_start));
                        }
                    }
                }
                reporter()->test_info->TearDownTestCase();
                has_failures = has_failures || reporter()->test_info->HasFailure();
                reporter()->EndCase(RunTime(case_start));
            }
            reporter()->End(RunTime(total_start));
            if (has_failures) return EXIT_FAILURE;
            return EXIT_SUCCESS;
        }
        
    };
    
    
    template <class C>
    class Factory : public Info {
        Test* CreateFixture() {
            return new C;
        }
        void SetUpTestCase() {
            C::SetUpTestCase();
        }
        void TearDownTestCase() {
            C::TearDownTestCase();
        }
    public:
        Factory(const char* test_case_name, const char* test_name, const char* type_name) :
        Info(test_case_name, test_name, type_name) {
            Runner::AddTest(this);
        }
    };

    
}
