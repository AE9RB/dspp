// benchtest - A benchmarking and unit testing framework.
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

namespace testing {

    template <typename T>
    void PrintTo(const T& val, ::std::ostream* os) {
        *os << val;
    }

    inline void PrintTo(const bool& val, ::std::ostream* os) {
        if (val) *os << "true";
        else *os << "false";
    }

    inline void PrintTo(const float& val, ::std::ostream* os) {
        *os << std::setprecision(::std::numeric_limits<float>::digits10) << val;
    }

    inline void PrintTo(const double& val, ::std::ostream* os) {
        *os << std::setprecision(::std::numeric_limits<double>::digits10) << val;
    }

    inline void PrintTo(const long double& val, ::std::ostream* os) {
        *os << std::setprecision(::std::numeric_limits<long double>::digits10) << val;
    }


    template <typename T>
    ::std::string PrintToString(const T& val) {
        ::std::stringstream ss;
        PrintTo(val, &ss);
        return ss.str();
    }

}
