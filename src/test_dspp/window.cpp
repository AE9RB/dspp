// dspp - Digital Signal Processing library for C++
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

#include <benchtest.hpp>
#include <dspp/window.hpp>

using namespace dspp;
using namespace std;

// These window tests do not validate the math. The only way to do
// that is to check against constants gathered from another program.
// Instead we check for symmetry on both odd and even sized functions.
// Also, we verify windows converge to zero at either their ends or one
// past the ends, depending on the type. This seems to catch all the
// human mistakes quite well.

template<typename T>
class Window : public testing::Test {
protected:

    // convenience for development
    void print(vector<T> &data) {
        cout << "[";
        for (auto i : data) cout << i << ", ";
        cout << "]\n";
    }
    void print(vector<T> &&data) {
        cout << "[";
        for (auto i : data) cout << i << ", ";
        cout << "]\n";
    }
    void print(function<vector<T>(size_t, bool)> fn) {
        cout << "symm8 ";
        print(fn(8, true));
        cout << "symm9 ";
        print(fn(9, true));
        cout << "peri8 ";
        print(fn(8, false));
        cout << "peri9 ";
        print(fn(9, false));
    }

    // Check middle and ends for each combination of
    // even/odd and symmetric/periodic
    void check(function<vector<T>(size_t, bool)> fn, T s8, T s9, T p8, T p9) {
        size_t last, mid;

        auto symm8 = fn(8, true);
        last = symm8.size()-1;
        mid = last/2;
        EXPECT_EQ(symm8[mid], symm8[mid+1]);
        EXPECT_EQ(symm8[0], symm8[last]);
        EXPECT_EQ(s8, symm8[last]);

        auto symm9 = fn(9, true);
        last = symm9.size()-1;
        mid = symm9.size()/2;
        EXPECT_EQ(symm9[mid-1], symm9[mid+1]);
        EXPECT_EQ(symm9[0], symm9[last]);
        EXPECT_EQ(s9, symm9[last]);

        auto peri8 = fn(8, false);
        last = peri8.size()-1;
        mid = peri8.size()/2;
        EXPECT_EQ(peri8[mid-1], peri8[mid+1]);
        EXPECT_EQ(peri8[1], peri8[last]);
        EXPECT_EQ(p8, peri8[last]);

        auto peri9 = fn(9, false);
        last = peri9.size()-1;
        mid = peri9.size()/2;
        EXPECT_EQ(peri9[mid-1], peri9[mid+1]);
        EXPECT_EQ(peri9[0], peri9[last]);
        EXPECT_EQ(p9, peri9[last]);

        // size of 1 should always have single value of 1
        EXPECT_EQ(1, fn(1, false)[0]);
        EXPECT_EQ(1, fn(1, true)[0]);
    }

    void rect() {
        check([](size_t size, bool symm) {
            return window::rect(vector<T>(size,1),symm);
        }, 1, 1, 1, 1);
    }

    void triang() {
        check([](size_t size, bool symm) {
            return window::triang(vector<T>(size,1),symm);
        }, 0.125, 0.2, 0.4, 0.2);
    }

    void bartlett() {
        check([](size_t size, bool symm) {
            return window::bartlett(vector<T>(size,1),symm);
        }, 0, 0, 0.25, 0);
    }

    void hann() {
        check([](size_t size, bool symm) {
            return window::hann(vector<T>(size,1),symm);
        }, 0, 0, 0.14644660940672627, 0);
    }

    void welch() {
        check([](size_t size, bool symm) {
            return window::welch(vector<T>(size,1),symm);
        }, 0.39506172839506171, 0.36, 0.64, 0.36);
    }

    void parzen() {
        check([](size_t size, bool symm) {
            return window::parzen(vector<T>(size,1),symm);
        }, 0.00390625, 0.0027434842249657101,
        0.074074074074074098, 0.0027434842249657101);
    }

    void bohman() {
        check([](size_t size, bool symm) {
            return window::bohman(vector<T>(size,1),symm);
        }, 0, 0, 0.048302383742639676, 0);
    }

    void chebyshev() {
        vector<T> w;
        w = window::chebyshev(vector<T>(1,1));
        EXPECT_EQ((T)1, w[0]);
        w = window::chebyshev(vector<T>(8,1));
        EXPECT_EQ((T)0.03638368090334488, w[7]);
        EXPECT_EQ((T)1, w[3]);
        EXPECT_EQ((T)1, w[4]);
        w = window::chebyshev(vector<T>(9,1));
        EXPECT_EQ((T)0.021827407475211173, w[8]);
        EXPECT_EQ((T)1, w[4]);
    }

};

TEST_T(Window, double, rect);
TEST_T(Window, double, hann);
TEST_T(Window, double, welch);
TEST_T(Window, double, parzen);
TEST_T(Window, double, triang);
TEST_T(Window, double, bartlett);
TEST_T(Window, double, bohman);
TEST_T(Window, double, chebyshev);

TEST_T(Window, float, rect);
TEST_T(Window, float, hann);
TEST_T(Window, float, welch);
TEST_T(Window, float, parzen);
TEST_T(Window, float, triang);
TEST_T(Window, float, bartlett);
TEST_T(Window, float, bohman);
TEST_T(Window, float, chebyshev);
