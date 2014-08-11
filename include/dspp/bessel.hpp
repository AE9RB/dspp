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

#ifndef DSPP_BESSEL_HPP
#define DSPP_BESSEL_HPP

#include <cmath>

namespace dspp {
namespace bessel {

template <typename T>
static T i0(T x)
{
    T y, ax = fabs(x);
    if ((ax=fabs(x)) < 3.75) {
        y = x/3.75;
        y = y*y;
        return 1.0+y*
               (3.5156229+y*
                (3.0899424+y*
                 (1.2067492+y*
                  (0.2659732+y*
                   (0.360768e-1+y*0.45813e-2)))));
    } else {
        y=3.75/ax;
        return (exp(ax)/sqrt(ax))*
               (0.39894228+y*
                (0.1328592e-1+y*
                 (0.225319e-2+y*
                  (-0.157565e-2+y*
                   (0.916281e-2+y*
                    (-0.2057706e-1+y*
                     (0.2635537e-1+y*
                      (-0.1647633e-1+y*0.392377e-2))))))));
    }
}

}
}

#endif