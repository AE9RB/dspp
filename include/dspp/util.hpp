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

#ifndef DSPP_UTIL_HPP
#define DSPP_UTIL_HPP

#include <complex>
#include <cmath>
#include <functional>
#include <type_traits>

namespace dspp {

// Mathematical constants use the same syntax as Boost library.
#define DSPP_DEFINE_MATH_CONSTANT(name, x)\
template<typename T=double> inline constexpr T name() {return x;}
/// @brief 3.141592653589793238462643383279502884e+00
DSPP_DEFINE_MATH_CONSTANT(pi, 3.141592653589793238462643383279502884e+00)
/// @brief 6.283185307179586476925286766559005768e+00
DSPP_DEFINE_MATH_CONSTANT(two_pi, 6.283185307179586476925286766559005768e+00)

/// @brief Fast fused multiply and accumulate. Returns x * y + z.
/// @details
/// A fused multiply and accumulate avoids the rounding that would normally
/// happen after multiplication. This offers greater precision.
/// Some microprocessors do not have an FMA instruction so \::std::fma()
/// will be significantly slower than (x*y+z).
/// Most DSP cases prefer speed over precision which this function ensures
/// by using (x*y+z) on processors without an FMA instruction.
/// If you require precision over speed, use \::std::fma() instead.
template<typename T1, typename T2, typename T3>
typename ::std::enable_if <
    ::std::is_same<double, decltype(T1()+T2()+T3())>::value,
    decltype(T1()+T2()+T3())
>::type
inline fmac(T1 x, T2 y, T3 z) {
    #ifdef FP_FAST_FMA
    return ::std::fma(x,y,z);
    #else
    return x*y+z;
    #endif
}

/// @brief Fast fused multiply and accumulate. Returns x * y + z.
template<typename T1, typename T2, typename T3>
typename ::std::enable_if <
    ::std::is_same<float, decltype(T1()+T2()+T3())>::value,
    decltype(T1()+T2()+T3())
>::type
inline fmac(T1 x, T2 y, T3 z) {
    #ifdef FP_FAST_FMAF
    return ::std::fma(x,y,z);
    #else
    return x*y+z;
    #endif
}

/// Function mapper to wrap filter and window algorithms.
template<typename T>
class Fmap : public std::iterator<std::input_iterator_tag, size_t>
{
private:
    Fmap(Fmap *b, size_t index)
        : _fn(b->_fn), _size(b->_size), _index(index) {}
    std::function<T(size_t)> _fn;
    size_t _size;
    size_t _index;
public:
    Fmap(size_t size, std::function<T(size_t)> &&fn)
        : _fn(fn), _size(size), _index(0) {}
    bool operator==(Fmap<T> const & rhs) const {
        return (_index == rhs._index);
    }
    bool operator!=(Fmap<T> const & rhs) const {
        return !operator==(rhs);
    }
    void operator++() {
        ++_index;
    }
    void operator--() {
        --_index;
    }
    T operator*() const {
        return _fn(_index);
    }
    T operator[](size_t const & x) const {
        return _fn(x);
    }
    size_t size() const {
        return _size;
    }
    size_t index() const {
        return _index;
    }
    Fmap begin() {
        return Fmap(this, 0);
    }
    Fmap end() {
        return Fmap(this, _size);
    }
};

} /* namespace dspp */


/// @brief Specialize complex multiplication.
/// @details <tt>\#include <dspp/util.hpp></tt>
///
/// The ISO specification for C++ requires that complex multiplication
/// check for NaN results and adjust for mathematical correctness.
/// This is a significant performance problem for DSP algorithms
/// which never use NaN or Inf. GCC has an option (<tt>-fcx-limited-range</tt>)
/// to disable this but clang currently does not. Since we are allowed
/// to specialize existing functions in std we can fix this for all
/// compilers. If parts of your application require the ISO rules then make
/// sure you are not including any dspp headers in that compilation unit.
namespace std {
#define DSPP_SPECIALIZE_COMPLEX_MULTIPLICATION(T1, T2) \
std::complex<decltype(T1()+T2())> \
inline operator*(const std::complex<T1>& z, const std::complex<T2>& w) \
{   return std::complex<decltype(T1()+T2())>( \
    z.real()*w.real() - z.imag()*w.imag(), \
    z.imag()*w.real() + z.real()*w.imag() \
);}
/// Limited range complex multiplication.
DSPP_SPECIALIZE_COMPLEX_MULTIPLICATION(float, float)
/// Limited range complex multiplication.
DSPP_SPECIALIZE_COMPLEX_MULTIPLICATION(double, float)
/// Limited range complex multiplication.
DSPP_SPECIALIZE_COMPLEX_MULTIPLICATION(float, double)
/// Limited range complex multiplication.
DSPP_SPECIALIZE_COMPLEX_MULTIPLICATION(double, double)
} /* namespace std */

#endif /* DSPP_UTIL_HPP */

