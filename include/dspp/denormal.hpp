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

#ifndef DSPP_DENORMAL_HPP
#define DSPP_DENORMAL_HPP

#ifdef __SSE__
#include <xmmintrin.h>
#endif

#ifdef __SSE3__
#include <pmmintrin.h>
#endif

namespace dspp {

/// @brief Disable special (slow) handling of denormal floating point types.
/// @details This is implemented uniquely for each CPU architecture by
/// setting processor flags. Generally, these flags are per-thread.
///
/// Intel x86-64: enables FTZ and DAZ flags
inline void fast_denormals(bool fast=true) {

    if (fast) {
        #ifdef __SSE__
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        #endif
        #ifdef __SSE3__
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
        #endif
    } else {
        #ifdef __SSE__
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_OFF);
        #endif
        #ifdef __SSE3__
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_OFF);
        #endif
    }

}

} /* namespace dspp */

#endif /* DSPP_DENORMAL_HPP */
