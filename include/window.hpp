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

#ifndef DSPP_WINDOW_HPP
#define DSPP_WINDOW_HPP

#include <cmath>
#include "util.hpp"

namespace dspp {
namespace window {

/// \f[
/// w(n) = 1.
/// \f]
template<typename T>
Fmap<T>
rect(size_t size, bool symm = true) {
    return Fmap<T>(size, [](size_t index)->T {
        return 1;
    });
}

/// \f[
/// w(n)=1 - \left|\frac{n-\frac{N-1}{2}}{\frac{N+1}{2}}\right|
/// \f]
/// Even symmetric windows:
/// \f[
/// w(n)=1 - \left|\frac{n-\frac{N-1}{2}}{\frac{N}{2}}\right|
/// \f]
template<typename T>
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
/// w(n) = \frac{2}{N-1} \left(
///        \frac{N-1}{2} - \left|n - \frac{N-1}{2}\right|
///        \right)
/// \f]
template<typename T>
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
template<typename T>
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
/// \f]
template<typename T>
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

template<typename T>
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
