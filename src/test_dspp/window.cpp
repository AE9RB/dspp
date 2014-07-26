// dspp - Digital signal processing for C++
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
#include "dspp/window.hpp"

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

    void check(function<Fmap<T>(size_t, bool)> fn, int zero) {
        size_t last, mid;

        auto evensymm = fn(16, true);
        last = evensymm.size()-1;
        mid = last/2;
        EXPECT_EQ(evensymm[mid], evensymm[mid+1]);
        EXPECT_EQ(evensymm[0], evensymm[last]);

        auto oddsymm = fn(17, true);
        last = oddsymm.size()-1;
        mid = oddsymm.size()/2;
        EXPECT_EQ(oddsymm[mid-1], oddsymm[mid+1]);
        EXPECT_EQ(oddsymm[0], oddsymm[last]);

        auto even = fn(16, false);
        last = even.size()-1;
        mid = even.size()/2;
        EXPECT_EQ(even[mid-1], even[mid+1]);
        EXPECT_EQ(even[1], even[last]);

        auto odd = fn(17, false);
        last = odd.size()-1;
        mid = odd.size()/2;
        EXPECT_EQ(odd[mid-1], odd[mid+1]);
        EXPECT_EQ(odd[0], odd[last]);

        // size of 1 should always have single value of 1
        EXPECT_EQ(1, fn(1, false)[0]);
        EXPECT_EQ(1, fn(1, true)[0]);

        if (zero == 0) {
            // converge to 0 at the ends
            EXPECT_EQ(0, evensymm[0]);
            EXPECT_EQ(0, oddsymm[0]);
            EXPECT_EQ(0, even[0]);
            EXPECT_EQ(0, odd[0]);
        } else if (zero == -1) {
            // converge to 0 one beyond the ends
            EXPECT_EQ(0, evensymm[evensymm.size()]);
            EXPECT_EQ(0, oddsymm[oddsymm.size()]);
            EXPECT_EQ(0, even[even.size()+1]);
            EXPECT_EQ(0, odd[odd.size()]);
        } // else don't check, special case
    }

    // convenience for development
    void print(Fmap<T> &&data) {
        std::cout << "[";
        for (auto i : data) std::cout << i << ", ";
        std::cout << "]\n";
    }
    void print(function<Fmap<T>(size_t, bool)> fn, size_t size) {
        std::cout << "symm=true  ";
        print(fn(size, true));
        std::cout << "symm=false ";
        print(fn(size, false));
    }

    void triang() {
        check([](size_t size, bool symm) {
            return window::triang<T>(size,symm);
        }, 1);

        // even symmetric never converges to zero
        auto evensymm = window::triang<T>(8, true);
        auto oddsymm = window::triang<T>(9, true);
        auto even = window::triang<T>(8, false);
        auto odd = window::triang<T>(9, false);
        EXPECT_EQ(0.125, evensymm[0]);
        EXPECT_EQ(0, oddsymm[oddsymm.size()]);
        EXPECT_EQ(0, even[even.size()+1]);
        EXPECT_EQ(0, odd[odd.size()]);
    }

    void bartlett() {
        check([](size_t size, bool symm) {
            return window::bartlett<T>(size,symm);
        }, 0);
    }

    void rect() {
        check([](size_t size, bool symm) {
            return window::rect<T>(size,symm);
        }, 1);
    }

    void hann() {
        check([](size_t size, bool symm) {
            return window::hann<T>(size,symm);
        }, 0);
    }

    void welch() {
        check([](size_t size, bool symm) {
            return window::welch<T>(size,symm);
        }, -1);
    }

    void parzen() {
        check([](size_t size, bool symm) {
            return window::parzen<T>(size,symm);
        }, 1);
        auto evensymm = window::parzen<T>(8, true);
        auto even = window::parzen<T>(8, false);
        EXPECT_EQ(0.00390625, evensymm[0]);
        EXPECT_EQ(0.00274348422496571, even[0]);
    }

};

TEST_T(Window, double, rect);
TEST_T(Window, double, hann);
TEST_T(Window, double, welch);
TEST_T(Window, double, parzen);
TEST_T(Window, double, triang);
TEST_T(Window, double, bartlett);
