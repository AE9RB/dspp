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

#ifndef DSPP_WINDOW_HPP
#define DSPP_WINDOW_HPP

#include <cmath>
#include "util.hpp"

namespace dspp {
/// @brief Window functions (also known as an apodization functions or
///        tapering functions).
///
/// A [window function](http://en.wikipedia.org/wiki/Window_function)
/// is a mathematical function that is zero-valued outside
/// of an interval. Typical applications include the design of
/// [finite impulse response (FIR) filters](http://en.wikipedia.org/wiki/Fir_filter)
/// and controlling [spectral leakage](http://en.wikipedia.org/wiki/Spectral_leakage)
/// in [Fourier analysis](http://en.wikipedia.org/wiki/Fourier_analysis).
/// When a discrete-time series is multiplied against the coefficients of
/// a window function everything outside the interval is zero-valued.
/// This gives the appearance of "looking through a window".
///
/// N represents the number of samples in a discrete-time window function
/// \f$w(n),\space 0 \le n \le N-1\f$. When N is an odd number, non-flat windows
/// have a singular maximum point. When N is even, they have a double maximum.
///
/// For spectral analysis, window functions are often required to be an even
/// length but still need a single maximum value. This is accomplished by
/// deleting the right-most coefficient. The dspp window functions can do this
/// automatically when the \p symm parameter is false. For example, a window of
/// \p size 1024 with \p symm false will be computed as if N is 1025.
///
/// <h2 class="groupheader"> Function Prototype</h2>
/// <tt>\#include <dspp/window.hpp></tt>
///
/// All window functions have the same prototype.
/// @tparam T Numeric type of generated points, e.g., float or double.
/// @retval Fmap<T> An object you can iterate over or randomly access.
/// @param size Number of samples to generate.
/// @param symm True generates a symmetric window for filter design.<br>
///             False generates a periodic window for spectral analysis.
namespace window {

/// \f[
/// w(n)=1
/// \f]
template<typename T = double>
Fmap<T>
rect(size_t size, bool symm = true) {
    return Fmap<T>(size, [](size_t index)->T {
        return 1;
    });
}

/// @headerfile "include/dspp/fft.hpp"
/// \f[
/// w(n)=1-\begin{cases}
/// \left|\frac{n-\frac{N-1}{2}}{\frac{N}{2}}\right|
/// & \text{if $N$ is even} \\
/// \left|\frac{n-\frac{N-1}{2}}{\frac{N+1}{2}}\right|
/// & \text{if $N$ is odd}
/// \end{cases}
/// \qquad 0 \leq n \leq N-1
/// \f]
template<typename T = double>
Fmap<T>
triang(size_t size, bool symm = true) {
    T len = size;
    bool odd = size & 1;
    if (!symm && !odd) ++len;
    T midm = (len-1) / 2;
    if (!symm || odd) ++len;
    T midp = len / 2;
    return Fmap<T>(size, [=](size_t index)->T {
        return 1 - fabs((index - midm) / midp);
    });
}

/// \f[
/// w(n) = 1- \left|\frac{n-\frac{N-1}{2}}{\frac{N-1}{2}}\right|
///        \qquad 0 \leq n \leq N-1
/// \f]
template<typename T = double>
Fmap<T>
bartlett(size_t size, bool symm = true) {
    if (size==1) return rect<T>(1);
    T len = size;
    bool odd = size & 1;
    if (!symm && !odd) ++len;
    T midm = (len - 1) / 2;
    return Fmap<T>(size, [=](size_t index)->T {
        return 1 - fabs((index - midm) / midm);
    });
}

/// \f[
/// w(n) = 0.5 - 0.5 \cos\left(\frac{2\pi{n}}{N-1}\right)
///        \qquad 0 \leq n \leq N-1
/// \f]
template<typename T = double>
Fmap<T>
hann(size_t size, bool symm = true) {
    if (size==1) return rect<T>(1);
    T len = size - 1;
    bool odd = size & 1;
    if (!symm && !odd) ++len;
    return Fmap<T>(size, [=](size_t index)->T {
        return .5 - .5 * cos((two_pi<T>()*index)/len);
    });
}

/// \f[
/// w(n)=1 - \left(\frac{n-\frac{N-1}{2}}{\frac{N+1}{2}}\right)^2
///          \qquad 0 \leq n \leq N-1
/// \f]
template<typename T = double>
Fmap<T>
welch(size_t size, bool symm = true) {
    T len = size;
    bool odd = size & 1;
    if (!symm && !odd) ++len;
    T midm = (len-1) / 2;
    T midp = (len+1) / 2;
    return Fmap<T>(size, [=](size_t index)->T {
        return 1 - pow(((index - midm)/midp), 2);
    });
}

/// \f[
/// w(n)=\begin{cases}
/// 1 - 6\left(\frac{|n|}{N/2}\right)^2 + 6\left(\frac{|n|}{N/2}\right)^3
/// & 0 \leq |n| \leq (N-1)/4
/// \\ 2\left(1-\frac{|n|}{N/2}\right)^3
/// & (N-1)/4 \lt |n| \le (N-1)/2
/// \end{cases}
/// \qquad -\frac{N-1}{2} \leq n \leq \frac{N-1}{2}
/// \f]
template<typename T = double>
Fmap<T>
parzen(size_t size, bool symm = true) {
    T len = size;
    bool odd = size & 1;
    if (!symm && !odd) ++len;
    T half = len / 2;
    T quad = half / 2;
    return Fmap<T>(size, [=](size_t index)->T {
        T i = fabs(index + 0.5 - half);
        if (i <= quad) {
            return 1 -
            6 * pow((i/half), 2) +
            6 * pow((i/half), 3);
        } else {
            return 2 * pow(1-(i/half), 3);
        }
    });
}

} /* namespace window */
} /* namespace dspp */

#endif /* DSPP_WINDOW_HPP */
