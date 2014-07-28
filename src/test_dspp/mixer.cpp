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
#include <dspp/mixer.hpp>

using namespace std;

template<typename T>
class Mixer : public testing::Test {
protected:

    // Starting with a signal at 0Hz, mix then count the
    // number of direction changes and zero crossings.
    void
    check(T rate, T freq) {
        vector<complex<T>> data(rate, complex<T>(0,1));
        dspp::Mixer<decltype(data)> mixer(rate, freq);
        complex<T> prev {data[0]};
        int dir {0};
        int dir_changes {0};
        int zero_crossings {0};
        mixer(data);
        for (auto &d : data) {
            if (d.real() < prev.real() && dir != -1) {
                dir = -1;
                ++dir_changes;
            } else if (d.real() > prev.real() && dir != 1) {
                dir = 1;
                ++dir_changes;
            }

            if (d.imag() > 0 && prev.imag() <= 0) ++zero_crossings;
            if (d.imag() < 0 && prev.imag() >= 0) ++zero_crossings;

            prev = d;
        }
        EXPECT_NEAR(freq*2, dir_changes, 1) << "For rate=" << rate << " and freq=" << freq;
        EXPECT_NEAR(freq*2, zero_crossings, 1) << "For rate=" << rate << " and freq=" << freq;
    }

    void
    correctness() {
        check(96000, 0);
        check(96000, 8000);
        check(8000, 1);
        check(8000, 4000);
    }

    void
    performance() {
        vector<complex<T>> data(96000, complex<T>(0.5,0.5));
        dspp::Mixer<decltype(data)> mixer;
        mixer.freq(1000);
        while (Benchmark()) {
            mixer(data);
        }
    }

};

TEST_T(Mixer, double, correctness)
TEST_T(Mixer, double, performance)
TEST_T(Mixer, float, correctness)
TEST_T(Mixer, float, performance)
