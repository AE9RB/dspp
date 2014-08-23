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

#include <dspp.hpp>
#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;
using namespace dspp;

#if defined(_WIN32)
#include <io.h>
#define POPEN _popen
#define PCLOSE _pclose
#define ACCESS _access
#else
#include <unistd.h>
#define POPEN popen
#define PCLOSE pclose
#define ACCESS access
#endif

class Gnuplot {
    FILE *gnuplot;
public:
    Gnuplot(const char* cmdname = GNUPLOT_CMD) {
        signal(SIGPIPE, SIG_IGN);
        gnuplot = POPEN(cmdname, "w");
        if (!gnuplot) exit(1);
    }

    ~Gnuplot() {
        if (gnuplot && PCLOSE(gnuplot)) exit(1);
    }

    template <typename T>
    Gnuplot& operator<<(const T& value) {
        stringstream ss;
        ss << value;
        fputs(ss.str().c_str(), gnuplot);
        return *this;
    }

    Gnuplot& operator<<(::std::ostream& (*basic_manipulator)(::std::ostream& stream)) {
        ::std::ostream& (*xendl)(::std::ostream& stream) = ::std::endl;
        ::std::ostream& (*xflush)(::std::ostream& stream) = ::std::flush;
        stringstream ss;
        ss << basic_manipulator;
        fputs(ss.str().c_str(), gnuplot);
        if (basic_manipulator == xendl || basic_manipulator == xflush) {
            fflush(gnuplot);
        }
        return *this;
    }

    Gnuplot& cmd(const std::string &cmdstr) {
        (*this) << cmdstr << endl;
        return *this;
    }
};

void plot_window(function<vector<double>(size_t, bool)> fn, const char* name, const char* filename) {
    Gnuplot p;
    auto w = fn(1025, true);
    auto aw = fn(8192, false);
    auto* o = new array<double, 810>;
    auto* a = new array<complex<double>, 8192>;

    double min = 0;
    for (auto i : w) if (i < min) min = i;

    double d = 0;
    for (auto i : aw) d+=i;

    // Usually, you see the Fourier transform of a window generated by placing
    // the window at the start of a zeroed array then running the dft.
    // However, it will not accurately show the side lobe attenuation.
    // This version sends an actual signal through the window and FFT.
    // The signal frequency is "wiggled" to gather 10 samples per bin.
    for (int j=0; j < 10; ++j) {
        for (int i=0; i < 8192; ++i) {
            double f = two_pi() * (i/(2.0 + (j-5)/20480.0));
            (*a)[i] = complex<double>(aw[i] * cos(f), aw[i] * sin(f));
        }
        FFT::dft(*a);
        for (int i=0; i < 81; ++i) {
            (*o)[i*10+j] = abs((*a)[i+4096-40]) / d;
        }
    }

    p.cmd("set term pngcairo size 768,240 font 'Lucida Grande,9'");
    p << "set output 'images/window_" << filename << ".png'" << endl;
    p.cmd("set title offset 0,-0.8 font ',11'");
    p.cmd("set size 0.5,1.03");
    p.cmd("set xtics offset 0,0.1");
    p.cmd("set xlabel offset 0,0.6");

    p.cmd("set grid");

    p.cmd("set multiplot");
    p.cmd("set ylabel offset 2.0,0.0");
    p.cmd("set linetype 1 lc rgb '#8B4513' lw 1.1");

    p << "set title '" << name << "'" << endl;
    p << "set yrange [" << min << ":1.04]" << endl;
    p.cmd("set ytics autofreq 0,0.1");
    p.cmd("set ylabel 'amplitude'");
    p.cmd("set xrange [0:1024]");
    p.cmd("set xtics ('0' 0, '' 128, '' 256, '' 384, '' 512, '' 640, '' 768, '' 896, 'N-1' 1024)");
    p.cmd("set xlabel 'samples'");
    p.cmd("plot '-' with lines notitle");
    for (auto v : w) {
        p << v << "\n";
    }
    p.cmd("e");

    p.cmd("set origin 0.5,0.00");
    p.cmd("set title 'Fourier transform'");
    p.cmd("set yrange [-130:5]");
    p.cmd("set ytics autofreq -130,10");
    p.cmd("set ylabel 'decibels'");
    p.cmd("set xrange [0:801]");
    p.cmd("set xtics ('-40' 0, '-30' 100, '-20' 200, '-10' 300, '0' 400, '10' 500, '20' 600, '30' 700, '40' 801)");
    p.cmd("set xlabel 'bins'");
    p.cmd("plot '-' with lines notitle");
    for (int i=5; i< 807; ++i) {
        p << 20*log10((*o)[i]) << "\n";
    }
    p.cmd("e");

    delete o;
    delete a;
}

int main() {
    if (ACCESS("images",0)) {
        cerr << "Directory \"images\" not found." << endl;
        exit(1);
    }

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::rect(w,symm);
        return w;
    }, "Rectangle window", "rect");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::triang(w,symm);
        return w;
    }, "Triangle window", "triang");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::bartlett(w,symm);
        return w;
    }, "Bartlett window", "bartlett");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::hann(w,symm);
        return w;
    }, "Hann window", "hann");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::welch(w,symm);
        return w;
    }, "Welch window", "welch");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::parzen(w,symm);
        return w;
    }, "Parzen window", "parzen");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::bohman(w,symm);
        return w;
    }, "Bohman window", "bohman");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::chebyshev(w);
        return w;
    }, "Chebyshev window (100dB)", "chebyshev");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::blackman(w,symm);
        return w;
    }, "Blackman window", "blackman");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::nuttall(w,symm);
        return w;
    }, "Nuttall window", "nuttall");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::blackmannuttall(w,symm);
        return w;
    }, "Blackman-Nuttall window", "blackmannuttall");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::blackmanharris(w,symm);
        return w;
    }, "Blackman-Harris window", "blackmanharris");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::flattop(w,symm);
        return w;
    }, "Flat top window", "flattop");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::barthann(w,symm);
        return w;
    }, "Bartlett–Hann window", "barthann");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::hamming(w,symm);
        return w;
    }, "Hamming window", "hamming");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::kaiser(w, 10, symm);
        return w;
    }, "Kaiser window (beta=10)", "kaiser");

    plot_window([](size_t size, bool symm) {
        vector<double> w(size,1);
        window::gaussian(w, 2.5, symm);
        return w;
    }, "Gaussian window (a=2.5)", "gaussian");

}
