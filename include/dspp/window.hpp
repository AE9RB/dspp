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
#include <vector>
#include "util.hpp"
#include "fft.hpp"
#include "bessel.hpp"

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
/// \p size() 1024 with \p symm false will be computed as if N is 1025.
///
/// <h2 class="groupheader"> Function Prototype</h2>
/// <tt>\#include <dspp/window.hpp></tt>
///
/// Window functions have the same prototype.
/// @retval T& A reference to the container that was windowed.
/// @param w A container of data for application of the window. To get the
///          window coefficients supply a container full of 1s. The container
///          must provide a forward iterator, T::value_type, and size().
/// @param symm True generates a symmetric window for filter design.<br>
///             False generates a periodic window for spectral analysis.
///
/// Some windows are generated with sine or FFT functions. If performance is
/// important then use a stored copy of the coefficients instead of generating
/// the window every time.
/// ~~~
///     // Obtaining window coefficients
///     auto w = dspp::window::hann(std::vector<float>(32,1));
/// ~~~
/// For creating FIR filters you can apply the window to an existing container.
/// ~~~
///     // Applying window to an existing container
///     std::array<double, 128> filter = some_fir_filter();
///     dspp::window::hann(filter);
/// ~~~
namespace window {

/// \f[
/// w(n)=1 \qquad 0 \leq n \leq N-1
/// \f]
/// @image html window_rect.png
template <class T>
T& rect(T&& w, bool symm=true) {
    (void)symm;
    return w;
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
/// @image html window_triang.png
template <class T>
T& triang(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    Tv len = w.size();
    bool odd = w.size() & 1;
    if (!symm && !odd) ++len;
    Tv midm = (len-1) / 2;
    if (!symm || odd) ++len;
    Tv midp = len / 2;
    Tv n = 0;
    for (auto &v : w) {
        v *= 1 - fabs((n - midm) / midp);
        ++n;
    }
    return w;
}

/// \f[
/// w(n) = 1- \left|\frac{n-\frac{N-1}{2}}{\frac{N-1}{2}}\right|
///        \qquad 0 \leq n \leq N-1
/// \f]
/// @image html window_bartlett.png
template <class T>
T& bartlett(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv midm = (len-1) / 2;
        Tv n = 0;
        for (auto &v : w) {
            v *= 1 - fabs((n - midm) / midm);
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n) = 0.5 - 0.5 \cos\left(\frac{2\pi{n}}{N-1}\right)
///        \qquad 0 \leq n \leq N-1
/// \f]
/// @image html window_hann.png
template <class T>
T& hann(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size() - 1;
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            v *= .5 - .5 * cos((two_pi<Tv>()*n)/len);
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n)=1 - \left(\frac{n-\frac{N-1}{2}}{\frac{N+1}{2}}\right)^2
///          \qquad 0 \leq n \leq N-1
/// \f]
/// @image html window_welch.png
template <class T>
T& welch(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    Tv len = w.size();
    bool odd = w.size() & 1;
    if (!symm && !odd) ++len;
    Tv midm = (len-1) / 2;
    Tv midp = (len+1) / 2;
    Tv n = 0;
    for (auto &v : w) {
        v *= 1 - pow(((n - midm)/midp), 2);
        ++n;
    }
    return w;
}

/// \f[
/// w(n)=\begin{cases}
/// 1 - 6\left(\frac{|n|}{N/2}\right)^2 + 6\left(\frac{|n|}{N/2}\right)^3
/// & 0 \leq |n| \leq N/4
/// \\ 2\left(1-\frac{|n|}{N/2}\right)^3
/// & N/4 \lt |n| \le N/2
/// \end{cases}
/// \qquad -\frac{N-1}{2} \leq n \leq \frac{N-1}{2}
/// \f]
/// @image html window_parzen.png
template <class T>
T& parzen(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    Tv len = w.size();
    bool odd = w.size() & 1;
    if (!symm && !odd) ++len;
    Tv half = len / 2;
    Tv quad = half / 2;
    Tv n = 0;
    for (auto &v : w) {
        Tv i = fabs(n + 0.5 - half);
        if (i <= quad) {
            v *= 1 -
                 6 * pow((i/half), 2) +
                 6 * pow((i/half), 3);
        } else {
            v *= 2 * pow(1-(i/half), 3);
        }
        ++n;
    }
    return w;
}

/// \f[
/// w(n) = (1-|n|)cos(\pi\left|n\right|)+\frac{1}{\pi}sin(\pi\left|n\right|)
/// \qquad -1 \leq n \leq 1
/// \f]
/// @image html window_bohman.png
template <class T>
T& bohman(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size() - 1;
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv half = len / 2;
        Tv n = 0;
        for (auto &v : w) {
            Tv x = fabs(n/half - 1);
            v *= (1 - x) * cos(pi<Tv>() * x) + 1.0 / pi<Tv>() * sin(pi<Tv>() * x);
            ++n;
        }
    }
    return w;
}

/// \f[
/// W(n) = \frac
/// {\cos\{N \cos^{-1}[\beta \cos(\frac{\pi k}{N})]\}}
/// {\cosh[N \cosh^{-1}(\beta)]}
/// \qquad
/// \beta = \cosh \left [\frac{1}{N} \cosh^{-1}(10^\frac{a}{20}) \right ]
/// \qquad -1 \leq n \leq 1
/// \f]
/// @image html window_chebyshev.png
template <class T>
T& chebyshev(T&& w, typename T::value_type a = 100) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        size_t len = w.size();
        bool odd = w.size() & 1;
        size_t order = len - 1.0;
        Tv beta = cosh(1.0 / order * acosh(pow(10, (abs(a) / 20))));
        ::std::vector<::std::complex<Tv>> k(len);
        for (size_t i = 0; i < len; ++i) {
            Tv x = beta * cos(pi<Tv>() * i / len);
            if (x>1) x = cosh(order * acosh(x));
            else if (x<-1) x = (1.0 - 2.0 * (order % 2)) * cosh(order * acosh(-x));
            else x =  cos(order * acos(x));
            if (odd)
                k[i] = x;
            else
                k[i] = x * exp(::std::complex<Tv>(0, pi<Tv>() / len * i));
        }
        if (!(len & (len-1))) FFT::dft(k);
        else FFT::czt(k.size(), k.data());
        size_t n = len / 2 + 1;
        Tv d;
        if (odd) d = k[0].real();
        else d = k[1].real();
        for (auto &v : w) {
            v *= k[n].real() / d;
            if (++n >= len) n = 0;
        }
    }
    return w;
}

/// \f[
/// w(n) = 0.42 - 0.5 \cos\left(\frac{2\pi n}{N-1}\right) +
/// 0.08 \cos\left(\frac{4\pi n}{N-1}\right)
/// \qquad 0 \leq n \leq N-1
/// \f]
/// @image html window_blackman.png
template <class T>
T& blackman(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            v *= 0.42 - 0.5 * cos(2.0 * pi<Tv>() * n / (len - 1)) +
                 0.08 * cos(4.0 * pi<Tv>() * n / (len - 1));
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n)=a_0 - a_1 \cos \left ( \frac{2 \pi n}{N-1} \right)+ a_2 \cos \left
/// ( \frac{4 \pi n}{N-1} \right)- a_3 \cos \left ( \frac{6 \pi n}{N-1} \right)
/// \qquad 0 \leq n \leq N-1
/// \\ a_0=0.355768;\quad a_1=0.487396;\quad a_2=0.144232;\quad a_3=0.012604
/// \f]
/// @image html window_nuttall.png
template <class T>
T& nuttall(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            v *= 0.355768 - 0.487396 * cos(2.0 * pi<Tv>() * n / (len - 1)) +
                 0.144232 * cos(4.0 * pi<Tv>() * n / (len - 1)) -
                 0.012604 * cos(6.0 * pi<Tv>() * n / (len - 1));
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n)=a_0 - a_1 \cos \left ( \frac{2 \pi n}{N-1} \right)+ a_2 \cos \left
/// ( \frac{4 \pi n}{N-1} \right)- a_3 \cos \left ( \frac{6 \pi n}{N-1} \right)
/// \qquad 0 \leq n \leq N-1
/// \\ a_0=0.3635819;\quad a_1=0.4891775;\quad a_2=0.1365995;\quad a_3=0.0106411
/// \f]
/// @image html window_blackmannuttall.png
template <class T>
T& blackmannuttall(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            v *= 0.3635819 - 0.4891775 * cos(2.0 * pi<Tv>() * n / (len - 1)) +
                 0.1365995 * cos(4.0 * pi<Tv>() * n / (len - 1)) -
                 0.0106411 * cos(6.0 * pi<Tv>() * n / (len - 1));
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n)=a_0 - a_1 \cos \left ( \frac{2 \pi n}{N-1} \right)+ a_2 \cos \left
/// ( \frac{4 \pi n}{N-1} \right)- a_3 \cos \left ( \frac{6 \pi n}{N-1} \right)
/// \qquad 0 \leq n \leq N-1
/// \\ a_0=0.35875;\quad a_1=0.48829;\quad a_2=0.14128;\quad a_3=0.01168
/// \f]
/// @image html window_blackmanharris.png
template <class T>
T& blackmanharris(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            v *= 0.35875 - 0.48829 * cos(2.0 * pi<Tv>() * n / (len - 1)) +
                 0.14128 * cos(4.0 * pi<Tv>() * n / (len - 1)) -
                 0.01168 * cos(6.0 * pi<Tv>() * n / (len - 1));
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n)=a_0 - a_1 \cos \left ( \frac{2 \pi n}{N-1} \right)+ a_2 \cos \left
/// ( \frac{4 \pi n}{N-1} \right)- a_3 \cos \left ( \frac{6 \pi n}{N-1} \right)
/// + a_4 \cos \left ( \frac{8 \pi n}{N-1} \right)
/// \qquad 0 \leq n \leq N-1
/// \\ a_0=0.21557895;\quad a_1=0.41663158;\quad a_2=0.277263158;\quad a_3=0.083578947;\quad a_4=0.006947368
/// \f]
/// @image html window_flattop.png
template <class T>
T& flattop(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            v *= 0.21557895 - 0.41663158 * cos(2.0 * pi<Tv>() * n / (len - 1)) +
                 0.277263158 * cos(4.0 * pi<Tv>() * n / (len - 1)) -
                 0.083578947 * cos(6.0 * pi<Tv>() * n / (len - 1)) +
                 0.006947368 * cos(8.0 * pi<Tv>() * n / (len - 1));
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n)=0.62 - 0.48 \left |\frac{n}{N-1}-\frac{1}{2} \right|
/// - 0.38 \cos \left (\frac{2 \pi n}{N-1}\right )
/// \qquad 0 \leq n \leq N-1
/// \f]
/// @image html window_barthann.png
template <class T>
T& barthann(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            Tv fac = fabs(n / (len - 1.0) - 0.5);
            v *= 0.62 - 0.48 * fac + 0.38 * cos(2 * pi<Tv>() * fac);
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n) = 0.54 - 0.46 \cos\left(\frac{2\pi{n}}{M-1}\right)
/// \qquad 0 \leq n \leq N-1
/// \f]
/// @image html window_hamming.png
template <class T>
T& hamming(T&& w, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            v *= 0.54 - 0.46 * cos(2.0 * pi<Tv>() * n / (len - 1));
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n) = \frac{I_0\left( \beta \sqrt{1-\frac{4n^2}{(N-1)^2}} \right)}
///        {I_0(\beta)}
/// \qquad 0 \leq n \leq N-1
/// \f]
/// @image html window_kaiser.png
template <class T>
T& kaiser(T&& w, typename T::value_type beta=10, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv alpha = (len - 1) / 2.0;
        Tv d = bessel::i0(beta);
        Tv n = 0;
        for (auto &v : w) {
            v *= bessel::i0(beta * sqrt(1 - pow(((n - alpha) / alpha), 2.0))) / d;
            ++n;
        }
    }
    return w;
}

/// \f[
/// w(n) = e^{ -\frac{1}{2}\left(\alpha\frac{n}{N/2}\right)^2 }
    /// \qquad -\frac{N-1}{2} \leq n \leq \frac{N-1}{2}
/// \f]
/// @image html window_gaussian.png
template <class T>
T& gaussian(T&& w, typename T::value_type a=2.5, bool symm=true) {
    typedef typename T::value_type Tv;
    if (w.size()>1) {
        Tv len = w.size();
        bool odd = w.size() & 1;
        if (!symm && !odd) ++len;
        Tv n = 0;
        for (auto &v : w) {
            Tv i = (n - (len - 1) / 2) * 2;
            v *= exp(-0.5 * pow((a / len * i), 2));
            ++n;
        }
    }
    return w;
}

} /* namespace window */
} /* namespace dspp */

#endif /* DSPP_WINDOW_HPP */
