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

#ifndef GNUPLOT_CMD
#define GNUPLOT_CMD "gnuplot"
#endif

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

template <typename F>
void plot_window(const char* name, const char* filename, F *fn) {
    Gnuplot p;
    auto w = fn(1025, true);
    auto aw = fn(8192, false);
    auto* o = new array<double, 810>;
    auto* a = new array<complex<double>, 8192>;

    double max = 0;
    for (auto i : aw) max+=i;

    for (int j=0; j < 10; ++j) {
        for (int i=0; i < 8192; ++i) {
            double f = two_pi() * (i/(2.0 + (j-5)/20480.0));
            (*a)[i] = complex<double>(aw[i] * cos(f), aw[i] * sin(f));
        }
        FFT::dft(*a);
        for (int i=0; i < 81; ++i) {
            (*o)[i*10+j] = abs((*a)[i+4096-40]) / max;
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
    p.cmd("set linetype 1 lc rgb '#8B4513' lw 1");

    p << "set title '" << name << " window'" << endl;
    p.cmd("set yrange [0.0:1.04]");
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
    p.cmd("set xrange [0:800]");
    p.cmd("set xtics ('-40' 0, '-30' 100, '-20' 200, '-10' 300, '0' 400, '10' 500, '20' 600, '30' 700, '40' 800)");
    p.cmd("set xlabel 'bins'");
    p.cmd("plot '-' with lines notitle");
    for (int i=5; i< 806; ++i) {
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

    plot_window("Rectangle", "rect", window::rect<double>);
    plot_window("Triangle", "triang", window::triang<double>);
    plot_window("Bartlett", "bartlett", window::bartlett<double>);
    plot_window("Hann", "hann", window::hann<double>);
    plot_window("Welch", "welch", window::welch<double>);
    plot_window("Parzen", "parzen", window::parzen<double>);
}
