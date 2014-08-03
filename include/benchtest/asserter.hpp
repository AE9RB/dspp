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

/// @cond BENCHTEST_IMPL
namespace testing {
    class AssertionResult {
    private:
        bool success;
        ::std::ostringstream* message;
    public:
        explicit AssertionResult(bool success) :
            success(success),
            message(new ::std::ostringstream) {
        }
        AssertionResult(const AssertionResult& other) :
            success(other.success),
            message(new ::std::ostringstream) {
            *message << other.str();
        }
        AssertionResult(AssertionResult&& other) :
            success(other.success),
            message(other.message) {
            other.message = nullptr;
        }

        ~AssertionResult() {
            delete message;
        }
        ::std::string str() const {
            return message->str();
        }
        template <typename T>
        AssertionResult& operator<<(const T& value) {
            PrintTo(value, message);
            return *this;
        }
        AssertionResult& operator<<(::std::ostream& (*basic_manipulator)(::std::ostream& stream)) {
            *message << basic_manipulator;
            return *this;
        }
        operator bool() {
            return success;
        }
    };
/// @endcond

    /// @ingroup benchtest_assert
    /// @brief Constructs a failure message that implicitly converts to bool.
    /// @details
    /// If you build a function that returns true or false to indicate
    /// passing or failing status, you can instead return AssertionFailure()
    /// or AssertionSuccess() with additional information.
    /// ~~~
    /// testing::AssertionResult IsEven(int n) {
    ///     if ((n % 2) == 0)
    ///         return ::testing::AssertionSuccess() << n << " is even";
    ///     else
    ///         return ::testing::AssertionFailure() << n << " is odd";
    /// }
    ///
    /// TEST(Numbers, IsEven){
    ///     EXPECT_TRUE(IsEven(4));
    ///     EXPECT_FALSE(IsEven(5));
    /// }
    /// ~~~
    inline AssertionResult AssertionFailure() {
        return AssertionResult(false);
    }

    /// @ingroup benchtest_assert
    /// @brief Opposite of AssertionFailure().
    inline AssertionResult AssertionSuccess() {
        return AssertionResult(true);
    }

/// @cond BENCHTEST_IMPL
    class ScopeTracer {
        struct TraceInfo {
            ::std::ostringstream message;
            const char* file;
            long line;
        };
        static ::std::vector<TraceInfo*>& trace() {
            static ::std::vector<TraceInfo*> trace;
            return trace;
        }
        TraceInfo* trace_info;
        bool do_pop = true;
    public:
        ScopeTracer(const char* file, long line) {
            trace().push_back(new TraceInfo());
            trace_info = trace().back();
            trace_info->file = file;
            trace_info->line = line;
        }
        ScopeTracer(ScopeTracer& other) :
            trace_info(other.trace_info) {
            other.do_pop = false;
        }
        ScopeTracer(ScopeTracer&& other) :
            trace_info(other.trace_info) {
            other.do_pop = false;
        }
        ~ScopeTracer() {
            if (do_pop) {
                delete trace().back();
                trace().pop_back();
            }
        }
        template <typename T>
        ScopeTracer& operator<<(const T& value) {
            PrintTo(value, &trace_info->message);
            return *this;
        }
        ScopeTracer& operator<<(::std::ostream& (*basic_manipulator)(::std::ostream& stream)) {
            trace_info->message << basic_manipulator;
            return *this;
        }
        
        static void Report() {
            if (!trace().empty()) {
                reporter()->Print("Scoped trace:");
                for (auto& t: trace()) {
                    reporter()->Trace(t->message.str(), t->file, t->line);
                }
            }
        }
    };

    class ResultReporter {
        AssertionResult* result;
        const char* file;
        long line;
    public:
        ResultReporter(AssertionResult& result, const char* file, int line, bool fatal) :
            result(&result),
            file(file),
            line(line) {
            if (fatal) ++reporter()->test_info->fatal_failure_count;
            else ++reporter()->test_info->nonfatal_failure_count;
        }
        const void operator+=(const AssertionResult& user) {
            auto user_str = user.str();

            auto result_str_empty = result->str().empty();
            if (!user_str.empty()) {
                if (!result_str_empty) *result << ::std::endl;
                *result << user_str;
            }
            else if (result_str_empty) {
                *result << "No message.";
            }
            ScopeTracer::Report();
            reporter()->Error(result->str(), file, line);
        }
    };

    inline AssertionResult ResultEq(const char* expected_expression,
                             const char* actual_expression,
                             const ::std::string& expected_value,
                             const ::std::string& actual_value) {
        auto msg = AssertionFailure();
        msg << "Value of: " << actual_expression;
        if (actual_value != actual_expression) {
            msg << ::std::endl << "  Actual: " << actual_value;
        }
        msg << ::std::endl << "Expected: " << expected_expression;
        if (expected_value != expected_expression) {
            msg << ::std::endl << "Which is: " << expected_value;
        }
        return msg;
    }

    template<typename T>
    ::std::string ResultPred(const char* expr, const T& value, bool& where) {
        ::std::ostringstream result;
        auto vs = PrintToString(value);
        if (vs != expr) {
            if (!where) {
                where = true;
                result << ", where";
            }
            result << ::std::endl << expr << " evaluates to " << vs;
        }
        return result.str();
    }

    namespace assert {

        template <typename Pred, typename T1>
        AssertionResult PRED(const char* pred_text, const char* e1,
                             Pred pred, const T1& v1) {
            if (pred(v1)) return AssertionSuccess();
            auto result = AssertionFailure();
            bool where = false;
            result << pred_text << "(" << e1 << ") evaluates to false";
            result << ResultPred(e1, v1, where);
            return result;
        }

        template <typename Pred, typename T1, typename T2>
        AssertionResult PRED(const char* pred_text, const char* e1, const char* e2,
                             Pred pred, const T1& v1, const T2& v2) {
            if (pred(v1, v2)) return AssertionSuccess();
            auto result = AssertionFailure();
            bool where = false;
            result << pred_text << "(" << e1 << ", " << e2 << ") evaluates to false";
            result << ResultPred(e1, v1, where);
            result << ResultPred(e2, v2, where);
            return result;
        }

        template <typename Pred, typename T1, typename T2, typename T3>
        AssertionResult PRED(const char* pred_text, const char* e1, const char* e2, const char* e3,
                             Pred pred, const T1& v1, const T2& v2, const T3& v3) {
            if (pred(v1, v2, v3)) return AssertionSuccess();
            auto result = AssertionFailure();
            bool where = false;
            result << pred_text << "(" << e1 << ", " << e2 << ", " << e3 << ") evaluates to false";
            result << ResultPred(e1, v1, where);
            result << ResultPred(e2, v2, where);
            result << ResultPred(e3, v3, where);
            return result;
        }

        template <typename T1, typename T2, typename T3>
        AssertionResult NEAR(const char* e1, const char* e2, const char* e3,
                             const T1& v1, const T2& v2, const T3& v3) {
            const auto diff = fabs(v1 - v2);
            if (diff <= v3) return AssertionSuccess();
            auto result = AssertionFailure();
            bool where = false;
            result << "The difference between " << e1 << " and " << e1
                   << " is " << diff << ", which exceeds " << v3;
            result << ResultPred(e1, v1, where);
            result << ResultPred(e2, v2, where);
            result << ResultPred(e3, v3, where);
            return result;
        }

        inline AssertionResult EQ(const char* e1, const char* e2, bool v1, AssertionResult v2) {
            if (v1 == v2) return AssertionSuccess();
            auto e_value = PrintToString(v1);
            auto a_value = ::std::string(v2 ? "true" : "false");
            auto message = v2.str();
            if (message.size()) {
                a_value += " (" + message + ")";
            }
            return ResultEq(e1, e2, e_value, a_value);
        }

        inline AssertionResult EQ(const char* e1, const char* e2, bool v1, bool v2) {
            if (v1 == v2) return AssertionSuccess();
            return ResultEq(e1, e2, PrintToString(v1), PrintToString(v2));
        }

        // Specialization for floating point equality.
        template<typename T>
        AssertionResult EQ(const char* e1, const char* e2, T v1, T v2,
                           typename ::std::enable_if<::std::is_floating_point<T>::value >::type* = 0) {
            auto limit = ::std::numeric_limits<T>::epsilon() * 4;
            if (fabs(v1 - v2) <= limit) return AssertionSuccess();
            ::std::stringstream vs1, vs2;
            vs1 << ::std::setprecision(::std::numeric_limits<T>::digits10 + 2) << v1;
            vs2 << ::std::setprecision(::std::numeric_limits<T>::digits10 + 2) << v2;
            return ResultEq(e1, e2, vs1.str(), vs2.str());
        }

        // Specialization for complex floating point equality.
        template<typename T>
        AssertionResult EQ(const char* e1, const char* e2, ::std::complex<T> v1, ::std::complex<T> v2,
                           typename ::std::enable_if<::std::is_floating_point<T>::value >::type* = 0) {
            auto limit = ::std::numeric_limits<T>::epsilon() * 4;
            if (fabs(v1.real() - v2.real()) <= limit && fabs(v1.imag() - v2.imag()) <= limit)
                return AssertionSuccess();
            ::std::stringstream vs1, vs2;
            vs1 << ::std::setprecision(::std::numeric_limits<T>::digits10 + 2) << v1;
            vs2 << ::std::setprecision(::std::numeric_limits<T>::digits10 + 2) << v2;
            return ResultEq(e1, e2, vs1.str(), vs2.str());
        }

#define BENCHTEST_COMPARE_PRED_(op_name, op) \
        template <typename T1, typename T2> \
        AssertionResult op_name(const char* e1, const char* e2, const T1& v1, const T2& v2) { \
            if (v1 op v2) return AssertionSuccess(); \
            auto result = AssertionFailure(); \
            result << "Expected: (" << e1 << ") " #op " (" << e2 << ")" \
                   << ::std::endl << "  Actual: " << PrintToString(v1) \
                   << " vs " << PrintToString(v2); \
            return result; \
        }
        BENCHTEST_COMPARE_PRED_(EQ, ==);
        BENCHTEST_COMPARE_PRED_(NE, !=);
        BENCHTEST_COMPARE_PRED_(LE, <=);
        BENCHTEST_COMPARE_PRED_(LT, <);
        BENCHTEST_COMPARE_PRED_(GE, >=);
        BENCHTEST_COMPARE_PRED_(GT, >);

    }
}
/// @endcond

#define BENCHTEST_AMBIGUOUS_ELSE_BLOCKER_ switch (0) case 0: default:

#define BENCHTEST_ASSERT_ \
BENCHTEST_AMBIGUOUS_ELSE_BLOCKER_ if (auto bt_result_ =

#define BENCHTEST_RESULT_REPORTER_(file, line, fatal) \
::testing::ResultReporter(bt_result_, file, line, fatal) += ::testing::AssertionFailure()

#define BENCHTEST_RESULT_NONFATAL_ \
);else BENCHTEST_RESULT_REPORTER_(__FILE__, __LINE__, false)

#define BENCHTEST_RESULT_FATAL_ \
);else return BENCHTEST_RESULT_REPORTER_(__FILE__, __LINE__, true)

#define BENCHTEST_TOKEN_(a, b) BENCHTEST_TOKEN_CONCAT_(a, b)
#define BENCHTEST_TOKEN_CONCAT_(a, b) a ## b

/// @ingroup benchtest_assert
/// @hideinitializer @def FAIL
/// @brief Unconditionally fail the current test fatally.
#define FAIL() \
BENCHTEST_ASSERT_ \
::testing::AssertionFailure() \
BENCHTEST_RESULT_FATAL_

/// @ingroup benchtest_assert
/// @hideinitializer @def ADD_FAILURE
/// @brief Unconditionally add a failure.
#define ADD_FAILURE() \
BENCHTEST_ASSERT_ \
::testing::AssertionFailure() \
BENCHTEST_RESULT_NONFATAL_

/// @ingroup benchtest_assert
/// @hideinitializer @def ADD_FAILURE_AT
/// @brief Unconditionally add a failure with custom filename and line number.
#define ADD_FAILURE_AT(file, line) \
BENCHTEST_ASSERT_ \
::testing::AssertionFailure() \
);else BENCHTEST_RESULT_REPORTER_(file, line, false)

/// @ingroup benchtest_assert
/// @hideinitializer @def SCOPED_TRACE
/// @brief Add scope information to failure reports.
/// @details Your tests may be implemented in a way that makes it
/// difficult to determine where a failure has occurred. For example,
/// you may be calling a subroutine multiple times in the same test.
/// Using SCOPED_TRACE is a fast way to add helpful information.
/// ~~~
/// TEST(ThisThat, OtherThing) {
///     {
///         SCOPED_TRACE() << "Check 1";  // The message "Check 1" will appear
///         check(1);                     // on failures in this scope.
///     }
///     // Leaving scope, no more extra message.
///     check(2);
/// }
/// ~~~
#define SCOPED_TRACE() \
auto BENCHTEST_TOKEN_(benchtest_tracer_, __LINE__) = ::testing::ScopeTracer(__FILE__, __LINE__)

#define EXPECT_PRED1(pred, v1) \
BENCHTEST_ASSERT_ \
::testing::assert::PRED(#pred, #v1, pred, v1) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_PRED1(pred, v1) \
BENCHTEST_ASSERT_ \
::testing::assert::PRED(#pred, #v1, pred, v1) \
BENCHTEST_RESULT_FATAL_
#define EXPECT_PRED_FORMAT1(pred_format, v1) \
BENCHTEST_ASSERT_ \
pred_format(#v1, v1) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_PRED_FORMAT1(pred_format, v1) \
BENCHTEST_ASSERT_ \
pred_format(#v1, v1) \
BENCHTEST_RESULT_FATAL_

#define EXPECT_PRED2(pred, v1, v2) \
BENCHTEST_ASSERT_ \
::testing::assert::PRED(#pred, #v1, #v2, pred, v1, v2) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_PRED2(pred, v1, v2) \
BENCHTEST_ASSERT_ \
::testing::assert::PRED(#pred, #v1, #v2, pred, v1, v2) \
BENCHTEST_RESULT_FATAL_
#define EXPECT_PRED_FORMAT2(pred_format, v1, v2) \
BENCHTEST_ASSERT_ \
pred_format(#v1, #v2, v1, v2) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_PRED_FORMAT2(pred_format, v1, v2) \
BENCHTEST_ASSERT_ \
pred_format(#v1, #v2, v1, v2) \
BENCHTEST_RESULT_FATAL_

#define EXPECT_PRED3(pred, v1, v2, v3) \
BENCHTEST_ASSERT_ \
::testing::assert::PRED(#pred, #v1, #v2, #v3, pred, v1, v2, v3) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_PRED3(pred, v1, v2, v3) \
BENCHTEST_ASSERT_ \
::testing::assert::PRED(#pred, #v1, #v2, #v3, pred, v1, v2, v3) \
BENCHTEST_RESULT_FATAL_
#define EXPECT_PRED_FORMAT3(pred_format, v1, v2, v3) \
BENCHTEST_ASSERT_ \
pred_format(#v1, #v2, #v3, v1, v2, v3) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_PRED_FORMAT3(pred_format, v1, v2, v3) \
BENCHTEST_ASSERT_ \
pred_format(#v1, #v2, #v3, v1, v2, v3) \
BENCHTEST_RESULT_FATAL_

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_TRUE
/// @brief Condition evaluates to true.
#define EXPECT_TRUE(condition) \
EXPECT_PRED_FORMAT2(::testing::assert::EQ, true, condition)
#define ASSERT_TRUE(condition) \
ASSERT_PRED_FORMAT2(::testing::assert::EQ, true, condition)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_FALSE
/// @brief Condition evaluates to false.
#define EXPECT_FALSE(condition) \
EXPECT_PRED_FORMAT2(::testing::assert::EQ, false, condition)
#define ASSERT_FALSE(condition) \
ASSERT_PRED_FORMAT2(::testing::assert::EQ, false, condition)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_NEAR
/// @brief Distance between val1 and val2 is less than or equal to abs_error.
#define EXPECT_NEAR(val1, val2, abs_error) \
EXPECT_PRED_FORMAT3(::testing::assert::NEAR, val1, val2, abs_error)
#define ASSERT_NEAR(val1, val2, abs_error) \
ASSERT_PRED_FORMAT3(::testing::assert::NEAR, val1, val2, abs_error)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_EQ
/// @brief Equality. (val1 == val2)
/// @details
/// A special case is triggered when both values are the same type
/// of floating point.  In this case the values are equal if they are
/// within 4 [ULPs] (http://en.wikipedia.org/wiki/Unit_in_the_last_place) (at 1.0).
/// Complex numbers of the same floating point type are also compared this way.
/// This is ideal for DSP values which typically range from -1...1. Consider using
/// \ref EXPECT_NEAR if you are working with significantly smaller or larger numbers.
#define EXPECT_EQ(val1, val2) \
EXPECT_PRED_FORMAT2(::testing::assert::EQ, val1, val2)
#define ASSERT_EQ(val1, val2) \
ASSERT_PRED_FORMAT2(::testing::assert::EQ, val1, val2)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_NE
/// @brief Inequality. (val1 != val2)
#define EXPECT_NE(val1, val2) \
EXPECT_PRED_FORMAT2(::testing::assert::NE, val1, val2)
#define ASSERT_NE(val1, val2) \
ASSERT_PRED_FORMAT2(::testing::assert::NE, val1, val2)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_LE
/// @brief Less than or equal. (val1 <= val2)
#define EXPECT_LE(val1, val2) \
EXPECT_PRED_FORMAT2(::testing::assert::LE, val1, val2)
#define ASSERT_LE(val1, val2) \
ASSERT_PRED_FORMAT2(::testing::assert::LE, val1, val2)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_LT
/// @brief Less than. (val1 < val2)
#define EXPECT_LT(val1, val2) \
EXPECT_PRED_FORMAT2(::testing::assert::LT, val1, val2)
#define ASSERT_LT(val1, val2) \
ASSERT_PRED_FORMAT2(::testing::assert::LT, val1, val2)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_GE
/// @brief Greater than or equal. (val1 >= val2)
#define EXPECT_GE(val1, val2) \
EXPECT_PRED_FORMAT2(::testing::assert::GE, val1, val2)
#define ASSERT_GE(val1, val2) \
ASSERT_PRED_FORMAT2(::testing::assert::GE, val1, val2)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_GT
/// @brief Greater than. (val1 > val2)
#define EXPECT_GT(val1, val2) \
EXPECT_PRED_FORMAT2(::testing::assert::GT, val1, val2)
#define ASSERT_GT(val1, val2) \
ASSERT_PRED_FORMAT2(::testing::assert::GT, val1, val2)

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_PRED1
/// @brief Use a predicate function which returns a \c bool.
/// @details Also available are <tt>EXPECT_PRED2(pred, v1, v2)</tt>
/// and <tt>EXPECT_PRED3(pred, v1, v2, v3)</tt>.
///
/// If you have a function that returns true or false,
/// or any value that can be implicitly converted to bool,
/// this will print its arguments for free.
/// ~~~
/// int HasCommonDenominatorGT1(int x, int y) {
///     if (y == 0) return x;
///     if (y == 1) return 0;
///     return HasCommonDenominatorGT1(y, x%y);
/// }
///
/// TEST(Numbers, Denominator){
///     EXPECT_PRED2(HasCommonDenominatorGT1, 100, 20);
/// }
/// ~~~

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_PRED_FORMAT1
/// @brief Use a predicate function which returns an \c AssertionResult.
/// @details Also available are <tt>EXPECT_PRED_FORMAT2(pred_format, v1, v2)</tt>
/// and <tt>EXPECT_PRED_FORMAT3(pred_format, v1, v2, v3)</tt>.
///
/// The most flexibility for generating failure messages comes from
/// defining predicate functions that return AssertionSuccess() or
/// AssertionFailure(). In fact, most of the \c EXPECT|ASSERT macros are
/// implemented this way.
///
/// For each variable, the predicate function is passed a text string of
/// the expression and the evaluated result. You may use value types or
/// reference types for the result parameters.
/// ~~~
/// int HasCommonDenominatorGT1(int x, int y) {
///     if (y == 0) return x;
///     if (y == 1) return 0;
///     return HasCommonDenominatorGT1(y, x%y);
/// }
///
/// testing::AssertionResult AssertCommonDenominator(
///         const char* x_expr, const char* y_expr, int x, int y) {
///    if (HasCommonDenominatorGT1(x,y)) return testing::AssertionSuccess();
///    return ::testing::AssertionFailure()
///        << x_expr << " and " << y_expr << " (" << x << " and " << y
///        << ") have no common divisor as determined by HasCommonDenominatorGT1(x,y)."
/// }
///
/// TEST(Numbers, Denominator){
///     EXPECT_PRED_FORMAT2(AssertCommonDenominator, 99, 20);
/// }
/// ~~~

#define BENCHTEST_NO_FATAL_FAILURE_(statement) \
BENCHTEST_ASSERT_ ::testing::AssertionSuccess()) { \
    auto bt_ffc_ = ::testing::reporter()->test_info->fatal_failure_count; \
    {statement;} \
    if (bt_ffc_ != ::testing::reporter()->test_info->fatal_failure_count) { \
        bt_result_ << "Expected: " #statement " doesn't generate new fatal failures." \
        << ::std::endl << "  Actual: it does."; \
        goto BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__); \
    } \
} else BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__): \
if (false

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_NO_FATAL_FAILURE
/// @brief Tests that \c statement doesn't add any new fatal failures.
/// @details When a function terminates on a fatal (<tt>ASSERT</tt>) failure
/// the caller will continue to run when it returns. You can
/// stop the test with \c ASSERT_NO_FATAL_FAILURE or use
/// \c EXPECT_NO_FATAL_FAILURE to provide more information.
/// ~~~
/// TEST(FooBar, SubRoutines) {
///     ASSERT_NO_FATAL_FAILURE(Foo());
///     EXPECT_NO_FATAL_FAILURE({
///         Bar();
///     }) << "More information";
/// }
/// ~~~
#define EXPECT_NO_FATAL_FAILURE(statement) \
BENCHTEST_NO_FATAL_FAILURE_(statement) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_NO_FATAL_FAILURE(statement) \
BENCHTEST_NO_FATAL_FAILURE_(statement) \
BENCHTEST_RESULT_FATAL_

#define BENCHTEST_CHECK_THROW_(statement, expected_exception) \
BENCHTEST_ASSERT_ ::testing::AssertionSuccess()) { \
    bool bt_caught_ = false; \
    try {statement;} \
    catch (expected_exception const&) {bt_caught_=true;} \
    catch (...) { \
        bt_result_ << "Expected: " #statement " throws an exception of type " \
        #expected_exception << ::std::endl << "  Actual: it throws a different type."; \
        goto BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__); \
    } \
    if (!bt_caught_) { \
        bt_result_ << "Expected: " #statement " throws an exception of type " \
        #expected_exception << ::std::endl << "  Actual: it throws nothing."; \
        goto BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__); \
    } \
} else BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__): \
if (false

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_THROW
/// @brief The \c statement must throw an exception of \c expected_exception type.
#define EXPECT_THROW(statement, expected_exception) \
BENCHTEST_CHECK_THROW_(statement, expected_exception) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_THROW(statement, expected_exception) \
BENCHTEST_CHECK_THROW_(statement, expected_exception) \
BENCHTEST_RESULT_FATAL_

#define BENCHTEST_CHECK_ANY_THROW_(statement) \
BENCHTEST_ASSERT_ ::testing::AssertionSuccess()) { \
    bool bt_caught_ = false; \
    try {statement;} \
    catch (...) {bt_caught_=true;} \
    if (!bt_caught_) { \
        bt_result_ << "Expected: " #statement " throws an exception." \
        << ::std::endl << "  Actual: it doesn't."; \
        goto BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__); \
    } \
} else BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__): \
if (false

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_ANY_THROW
/// @brief The \c statement must throw an exception of any type.
#define EXPECT_ANY_THROW(statement) \
BENCHTEST_CHECK_ANY_THROW_(statement) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_ANY_THROW(statement) \
BENCHTEST_CHECK_ANY_THROW_(statement) \
BENCHTEST_RESULT_FATAL_

#define BENCHTEST_CHECK_NO_THROW_(statement) \
BENCHTEST_ASSERT_ ::testing::AssertionSuccess()) { \
    bool bt_caught_ = false; \
    try {statement;} \
    catch (...) {bt_caught_=true;} \
    if (bt_caught_) { \
        bt_result_ << "Expected: " #statement " doesn't throw an exception." \
        << ::std::endl << "  Actual: it does."; \
        goto BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__); \
    } \
} else BENCHTEST_TOKEN_(benchtest_fatal_, __LINE__): \
if (false

/// @ingroup benchtest_assert
/// @hideinitializer @def EXPECT_NO_THROW
/// @brief The \c statement must not throw any exception.
#define EXPECT_NO_THROW(statement) \
BENCHTEST_CHECK_NO_THROW_(statement) \
BENCHTEST_RESULT_NONFATAL_
#define ASSERT_NO_THROW(statement) \
BENCHTEST_CHECK_NO_THROW_(statement) \
BENCHTEST_RESULT_FATAL_
