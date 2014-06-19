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

#ifndef dspp_util_hpp
#define dspp_util_hpp

#include <complex>
#include <cmath>

namespace dspp {

    // Mathematical constants use the same syntax as Boost library.

#define DSPP_DEFINE_MATH_CONSTANT(name, x)\
template<typename T> inline constexpr T name() {return x;}

    DSPP_DEFINE_MATH_CONSTANT(pi, 3.141592653589793238462643383279502884e+00)
    DSPP_DEFINE_MATH_CONSTANT(two_pi, 6.283185307179586476925286766559005768e+00)

    // Limited range (fast) multiplication of complex numbers.
    // Does not test for (nan,nan) results like the std library.
    template<class T1, class T2>
    std::complex<decltype(T1()+T2())>
    inline mul(const std::complex<T1>& z, const std::complex<T2>& w)
    {
        return std::complex<decltype(T1()+T2())>(
                   z.real()*w.real() - z.imag()*w.imag(),
                   z.imag()*w.real() + z.real()*w.imag()
               );
    }

    // Fast multiply and accumulate. Some platforms can perform this with a single
    // instruction that doesn't lose precision in any intermediate result.
    // Most DSP cases prefer speed over precision which this function ensures.
    template<typename T1, typename T2, typename T3>
    decltype(T1()+T2()+T3())
    inline fma(T1 x, T2 y, T3 z) {
        #ifdef FP_FAST_FMA
        return std::fma(x,y,z);
        #else
        return x*y+z;
        #endif
    }

}
#endif
