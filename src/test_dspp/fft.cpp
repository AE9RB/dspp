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
#include <dspp/fft.hpp>

const std::array<std::complex<double>, 16> ref0 {{
        std::complex<double>(-0.82993510256513270,0.78322255460971535),
        std::complex<double>(-0.62062045620071216,-0.20398322370742217),
        std::complex<double>(0.48702490306452950,0.12077985630401211),
        std::complex<double>(0.61913330685474266,0.02342510560093802),
        std::complex<double>(0.99016909661405061,0.93322272660158068),
        std::complex<double>(-0.14789834511540456,0.30599745382135302),
        std::complex<double>(0.92306621915157949,0.71597467817430172),
        std::complex<double>(-0.41194770159675098,-0.17071084234348244),
        std::complex<double>(0.02978581035346006,0.57956906405743114),
        std::complex<double>(0.08854560347058538,-0.81274017619000083),
        std::complex<double>(-0.13548094921372478,0.68985487733912110),
        std::complex<double>(0.54569292817085513,-0.61628209105191778),
        std::complex<double>(0.56073352395029885,-0.63731363839781507),
        std::complex<double>(0.15828299873972318,-0.37173711567411705),
        std::complex<double>(-0.17603078925720173,0.98461092905467384),
        std::complex<double>(-0.67215518569117150,-0.33030366956422297)
    }
};

const std::array<std::complex<double>, 16> ref1 {{
        std::complex<double>(1.40836586072972647,1.99358648863414878),
        std::complex<double>(0.54234607399180934,-0.86832671182654908),
        std::complex<double>(-5.95220210043290621,-0.35178684495262602),
        std::complex<double>(-3.20159619686836550,0.16091220573396947),
        std::complex<double>(-0.33641761978400742,-0.25010608715771332),
        std::complex<double>(2.56755504427141545,1.32271441382052957),
        std::complex<double>(-0.22468069674494817,0.07223943121245702),
        std::complex<double>(0.08659485141171297,1.62441635509585947),
        std::complex<double>(2.29029956346599217,6.34625560685189249),
        std::complex<double>(-0.97941948863275696,-2.44283615041583202),
        std::complex<double>(-0.52980347229082314,3.27653485796653410),
        std::complex<double>(-2.16961932435643590,-0.90456361049029954),
        std::complex<double>(-0.35923449100100391,-1.45493318084468060),
        std::complex<double>(0.71278017869274402,1.08532011997598121),
        std::complex<double>(-2.69752138163541089,1.27054267762715867),
        std::complex<double>(-4.43640844185886429,1.65159130252461450)
    }
};

const std::array<std::complex<double>, 8> ref2 {{
        std::complex<double>(1.00899192020690176,2.50792830906099606),
        std::complex<double>(-3.70198435453455810,-0.60666385918657428),
        std::complex<double>(-1.00055716131071581,1.85539515330709071),
        std::complex<double>(-0.38489594400606586,-0.48297521012410377),
        std::complex<double>(2.13165831232315206,2.59847132231822364),
        std::complex<double>(-1.12861368756438774,1.17874614737694361),
        std::complex<double>(-1.49915709502366634,-0.09601365984112620),
        std::complex<double>(-2.06492281061172145,-0.68910776603372681)
    }
};


template<typename T>
class FFTcorrectness : public testing::Test {
    std::array<std::complex<T>, 16> test16_in;
    std::array<std::complex<T>, 16> test16_out;
    std::array<std::complex<T>, 8> test8_in;
    std::array<std::complex<T>, 8> test8_out;
protected:
    void czt() {
        auto limit = ::std::numeric_limits<T>::epsilon() * 128;
        for (size_t i=0; i<16; ++i) {
            test16_in[i] = ref0[i];
            if (i<8) test8_in[i] = ref0[i];
        }
        dspp::FFT::czt(test16_in.size(), test16_in.data());
        for (size_t i=0; i<16; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_NEAR(ref1[i].real(), test16_in[i].real(), limit);
            ASSERT_NEAR(ref1[i].imag(), test16_in[i].imag(), limit);
        }
        dspp::FFT::czt(test8_in.size(), test8_in.data());
        for (size_t i=0; i<8; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_NEAR(ref2[i].real(), test8_in[i].real(), limit);
            ASSERT_NEAR(ref2[i].imag(), test8_in[i].imag(), limit);
        }
    }
    void dft() {
        for (size_t i=0; i<16; ++i) {
            test16_in[i] = ref0[i];
            if (i<8) test8_in[i] = ref0[i];
        }
        dspp::FFT::dft(test16_in, test16_out);
        for (size_t i=0; i<16; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_EQ(std::complex<T>(ref1[i]), test16_out[i]);
        }
        dspp::FFT::dft(test16_in);
        for (size_t i=0; i<16; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_EQ(std::complex<T>(ref1[i]), test16_in[i]);
        }
        dspp::FFT::dft(test8_in, test8_out);
        for (size_t i=0; i<8; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_EQ(std::complex<T>(ref2[i]), test8_out[i]);
        }
        dspp::FFT::dft(test8_in);
        for (size_t i=0; i<8; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_EQ(std::complex<T>(ref2[i]), test8_in[i]);
        }
    }
    void idft() {
        for (size_t i=0; i<16; ++i) {
            test16_in[i] = ref1[i];
            if (i<8) test8_in[i] = ref2[i];
        }
        dspp::FFT::dft(test16_in, test16_out, 1);
        for (size_t i=0; i<16; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_EQ(std::complex<T>(ref0[i]), test16_out[i] / std::complex<T>(16));
        }
        dspp::FFT::dft(test16_in, 1);
        for (size_t i=0; i<16; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_EQ(std::complex<T>(ref0[i]), test16_in[i] / std::complex<T>(16));
        }
        dspp::FFT::dft(test8_in, test8_out, 1);
        for (size_t i=0; i<8; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_EQ(std::complex<T>(ref0[i]), test8_out[i] / std::complex<T>(8));
        }
        dspp::FFT::dft(test8_in, 1);
        for (size_t i=0; i<8; ++i) {
            SCOPED_TRACE() << "i=" << i;
            ASSERT_EQ(std::complex<T>(ref0[i]), test8_in[i] / std::complex<T>(8));
        }
    }
};

TEST_T(FFTcorrectness, double, czt)
TEST_T(FFTcorrectness, float, czt)
TEST_T(FFTcorrectness, double, dft)
TEST_T(FFTcorrectness, float, dft)
TEST_T(FFTcorrectness, double, idft)
TEST_T(FFTcorrectness, float, idft)


template<typename T>
class FFTperformance : public testing::Test {
    std::array<std::complex<T>, 8192> *data;
    void SetUp() {
        for (size_t i=0; i<8192; ++i) {
            (*data)[i] = ref0[i%16];
        }
    }
protected:
    FFTperformance() :
        data(new std::array<std::complex<T>, 8192>) {
    }
    ~FFTperformance() {
        delete data;
    }
    void dft() {
        while (Benchmark()) {
            dspp::FFT::dft(*data);
        }
    }
};

TEST_T(FFTperformance, double, dft)
TEST_T(FFTperformance, float, dft)
