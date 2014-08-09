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

#ifndef DSPP_FFT_HPP
#define DSPP_FFT_HPP

#include <complex>
#include <array>
#include <vector>
#include <cmath>
#include <cassert>

namespace dspp {

/// @cond DSPP_IMPL

// FFT_Impl is based on fft4g.c
// Copyright (C) 1996-2001 Takuya OOURA
// http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
//
// Changes from the original:
//  * C++ template for the floating point type.
//  * Assert for power of two checking.
//  * All planning and work areas are fully automatic.
//  * Added out-of-place transforms for some functions.
template <typename T>
class FFT_Impl {
    static ::std::vector<::std::vector<size_t>> &it;
    static size_t wt_n;
    static ::std::vector<T> wt;
    static size_t ct_n;
    static ::std::vector<T> ct;
public:

    ///-------- Complex DFT (Discrete Fourier Transform) --------
    /// [definition]
    ///     <case1>
    ///         X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
    ///     <case2>
    ///         X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n
    ///     (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
    /// [usage]
    ///     <case1>
    ///         cdft(2*n, 1, a);
    ///     <case2>
    ///         cdft(2*n, -1, a);
    /// [parameters]
    ///     2*n            :data length (int)
    ///                     n >= 1, n = power of 2
    ///     a[0...2*n-1]   :input/output data (double *)
    ///                     input data
    ///                         a[2*j] = Re(x[j]),
    ///                         a[2*j+1] = Im(x[j]), 0<=j<n
    ///                     output data
    ///                         a[2*k] = Re(X[k]),
    ///                         a[2*k+1] = Im(X[k]), 0<=k<n
    /// [remark]
    ///     Inverse of
    ///         cdft(2*n, -1, a);
    ///     is
    ///         cdft(2*n, 1, a);
    ///         for (j = 0; j <= 2 * n - 1; j++) {
    ///             a[j] *= 1.0 / n;
    ///         }
    static void cdft(size_t n, int isgn, T *a)
    {
        powcheck(n);
        if (n > (wt_n << 2)) {
            makewt(n >> 2);
        }
        if (n > 4) {
            if (isgn >= 0) {
                bitrv2(n, a);
                cftfsub(n, a);
            } else {
                bitrv2conj(n, a);
                cftbsub(n, a);
            }
        } else if (n == 4) {
            cftfsub(n, a);
        }
    }

    static void cdft(size_t n, int isgn, const T *a, T *b)
    {
        powcheck(n);
        if (n > (wt_n << 2)) {
            makewt(n >> 2);
        }
        if (n > 4) {
            if (isgn >= 0) {
                bitrv2(n, a, b);
                cftfsub(n, b);
            } else {
                bitrv2conj(n, a, b);
                cftbsub(n, b);
            }
        } else {
            for (size_t i = 0; i < n; ++i) b[i] = a[i];
            if (n == 4) cftfsub(n, b);
        }
    }

    ///-------- Real DFT / Inverse of Real DFT --------
    /// [definition]
    ///     <case1> RDFT
    ///         R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
    ///         I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
    ///     <case2> IRDFT (excluding scale)
    ///         a[k] = (R[0] + R[n/2]*cos(pi*k))/2 +
    ///                sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) +
    ///                sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
    /// [usage]
    ///     <case1>
    ///         rdft(n, 1, a);
    ///     <case2>
    ///         rdft(n, -1, a);
    /// [parameters]
    ///     n              :data length (int)
    ///                     n >= 2, n = power of 2
    ///     a[0...n-1]     :input/output data (double *)
    ///                     <case1>
    ///                         output data
    ///                             a[2*k] = R[k], 0<=k<n/2
    ///                             a[2*k+1] = I[k], 0<k<n/2
    ///                             a[1] = R[n/2]
    ///                     <case2>
    ///                         input data
    ///                             a[2*j] = R[j], 0<=j<n/2
    ///                             a[2*j+1] = I[j], 0<j<n/2
    ///                             a[1] = R[n/2]
    /// [remark]
    ///     Inverse of
    ///         rdft(n, 1, a);
    ///     is
    ///         rdft(n, -1, a);
    ///         for (j = 0; j <= n - 1; j++) {
    ///             a[j] *= 2.0 / n;
    ///         }
    static void rdft(size_t n, int isgn, T *a)
    {
        size_t nw, nc;
        T xi;

        powcheck(n);
        nw = wt_n;
        if (n > (nw << 2)) {
            nw = n >> 2;
            makewt(nw);
        }
        nc = ct_n;
        if (n > (nc << 2)) {
            nc = n >> 2;
            makect(nc);
        }
        if (isgn >= 0) {
            if (n > 4) {
                bitrv2(n, a);
                cftfsub(n, a);
                rftfsub(n, a, nc);
            } else if (n == 4) {
                cftfsub(n, a);
            }
            xi = a[0] - a[1];
            a[0] += a[1];
            a[1] = xi;
        } else {
            a[1] = 0.5 * (a[0] - a[1]);
            a[0] -= a[1];
            if (n > 4) {
                rftbsub(n, a, nc);
                bitrv2(n, a);
                cftbsub(n, a);
            } else if (n == 4) {
                cftfsub(n, a);
            }
        }
    }

    static void rdft(size_t n, int isgn, const T *a, T *b)
    {
        size_t nw, nc;
        T xi;

        powcheck(n);
        nw = wt_n;
        if (n > (nw << 2)) {
            nw = n >> 2;
            makewt(nw);
        }
        nc = ct_n;
        if (n > (nc << 2)) {
            nc = n >> 2;
            makect(nc);
        }
        if (isgn >= 0) {
            if (n > 4) {
                bitrv2(n, a, b);
                cftfsub(n, b);
                rftfsub(n, b, nc);
            } else {
                for (size_t i = 0; i < n; ++i) b[i] = a[i];
                if (n == 4) cftfsub(n, b);
            }
            xi = b[0] - b[1];
            b[0] += b[1];
            b[1] = xi;
        } else {
            b[1] = 0.5 * (a[0] - a[1]);
            b[0] -= b[1];
            if (n > 4) {
                rftbsub(n, a, b, nc);
                bitrv2(n, b);
                cftbsub(n, b);
            } else {
                for (size_t i = 2; i < n; ++i) b[i] = a[i];
                if (n == 4) cftfsub(n, b);
            }
        }
    }

    ///-------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
    /// [definition]
    ///     <case1> IDCT (excluding scale)
    ///         C[k] = sum_j=0^n-1 a[j]*cos(pi*j*(k+1/2)/n), 0<=k<n
    ///     <case2> DCT
    ///         C[k] = sum_j=0^n-1 a[j]*cos(pi*(j+1/2)*k/n), 0<=k<n
    /// [usage]
    ///     <case1>
    ///         ddct(n, 1, a);
    ///     <case2>
    ///         ddct(n, -1, a);
    /// [parameters]
    ///     n              :data length (int)
    ///                     n >= 2, n = power of 2
    ///     a[0...n-1]     :input/output data (double *)
    ///                     output data
    ///                         a[k] = C[k], 0<=k<n
    /// [remark]
    ///     Inverse of
    ///         ddct(n, -1, a, ip, w);
    ///     is
    ///         a[0] *= 0.5;
    ///         ddct(n, 1, a, ip, w);
    ///         for (j = 0; j <= n - 1; j++) {
    ///             a[j] *= 2.0 / n;
    ///         }
    static void ddct(int n, int isgn, double *a)
    {
        int j, nw, nc;
        double xr;

        powcheck(n);
        nw = wt_n;
        if (n > (nw << 2)) {
            nw = n >> 2;
            makewt(nw);
        }
        nc = ct_n;
        if (n > nc) {
            nc = n;
            makect(nc);
        }
        if (isgn < 0) {
            xr = a[n - 1];
            for (j = n - 2; j >= 2; j -= 2) {
                a[j + 1] = a[j] - a[j - 1];
                a[j] += a[j - 1];
            }
            a[1] = a[0] - xr;
            a[0] += xr;
            if (n > 4) {
                rftbsub(n, a, nc);
                bitrv2(n, a);
                cftbsub(n, a);
            } else if (n == 4) {
                cftfsub(n, a);
            }
        }
        dctsub(n, a, nc);
        if (isgn >= 0) {
            if (n > 4) {
                bitrv2(n, a);
                cftfsub(n, a);
                rftfsub(n, a, nc);
            } else if (n == 4) {
                cftfsub(n, a);
            }
            xr = a[0] - a[1];
            a[0] += a[1];
            for (j = 2; j < n; j += 2) {
                a[j - 1] = a[j] - a[j + 1];
                a[j] += a[j + 1];
            }
            a[n - 1] = xr;
        }
    }

    ///-------- DST (Discrete Sine Transform) / Inverse of DST --------
    /// [definition]
    ///     <case1> IDST (excluding scale)
    ///         S[k] = sum_j=1^n A[j]*sin(pi*j*(k+1/2)/n), 0<=k<n
    ///     <case2> DST
    ///         S[k] = sum_j=0^n-1 a[j]*sin(pi*(j+1/2)*k/n), 0<k<=n
    /// [usage]
    ///     <case1>
    ///         ddst(n, 1, a, ip, w);
    ///     <case2>
    ///         ddst(n, -1, a, ip, w);
    /// [parameters]
    ///     n              :data length (int)
    ///                     n >= 2, n = power of 2
    ///     a[0...n-1]     :input/output data (double *)
    ///                     <case1>
    ///                         input data
    ///                             a[j] = A[j], 0<j<n
    ///                             a[0] = A[n]
    ///                         output data
    ///                             a[k] = S[k], 0<=k<n
    ///                     <case2>
    ///                         output data
    ///                             a[k] = S[k], 0<k<n
    ///                             a[0] = S[n]
    /// [remark]
    ///     Inverse of
    ///         ddst(n, -1, a, ip, w);
    ///     is
    ///         a[0] *= 0.5;
    ///         ddst(n, 1, a, ip, w);
    ///         for (j = 0; j <= n - 1; j++) {
    ///             a[j] *= 2.0 / n;
    ///         }
    static void ddst(int n, int isgn, double *a)
    {
        int j, nw, nc;
        double xr;

        powcheck(n);
        nw = wt_n;
        if (n > (nw << 2)) {
            nw = n >> 2;
            makewt(nw);
        }
        nc = ct_n;
        if (n > nc) {
            nc = n;
            makect(nc);
        }
        if (isgn < 0) {
            xr = a[n - 1];
            for (j = n - 2; j >= 2; j -= 2) {
                a[j + 1] = -a[j] - a[j - 1];
                a[j] -= a[j - 1];
            }
            a[1] = a[0] + xr;
            a[0] -= xr;
            if (n > 4) {
                rftbsub(n, a, nc);
                bitrv2(n, a);
                cftbsub(n, a);
            } else if (n == 4) {
                cftfsub(n, a);
            }
        }
        dstsub(n, a, nc);
        if (isgn >= 0) {
            if (n > 4) {
                bitrv2(n, a);
                cftfsub(n, a);
                rftfsub(n, a, nc);
            } else if (n == 4) {
                cftfsub(n, a);
            }
            xr = a[0] - a[1];
            a[0] += a[1];
            for (j = 2; j < n; j += 2) {
                a[j - 1] = -a[j] - a[j + 1];
                a[j] -= a[j + 1];
            }
            a[n - 1] = -xr;
        }
    }

private:

    static void powcheck(size_t n)
    {
        assert((n>=2) && !(n & (n-1))); // FFT size must be power of 2
    }

    static void makewt(size_t nw)
    {
        wt_n = nw;
        wt.reserve(nw);
        wt.resize(nw);
        wt.shrink_to_fit();
        size_t nwh = nw >> 1;
        if (nw > 2) {
            double delta = atan(1.0) / nwh;
            wt[0] = 1;
            wt[1] = 0;
            wt[nwh] = cos(delta * nwh);
            wt[nwh + 1] = wt[nwh];
            if (nwh > 2) {
                for (size_t j = 2; j < nwh; j += 2) {
                    double x = cos(delta * j);
                    double y = sin(delta * j);
                    wt[j] = x;
                    wt[j + 1] = y;
                    wt[nw - j] = y;
                    wt[nw - j + 1] = x;
                }
                bitrv2(nw, wt.data());
            }
        }
    }

    static void makect(size_t nc)
    {
        ct_n = nc;
        ct.reserve(nc);
        ct.resize(nc);
        ct.shrink_to_fit();
        size_t nch = nc >> 1;
        if (nc > 1) {
            double delta = atan(1.0) / nch;
            ct[0] = cos(delta * nch);
            ct[nch] = 0.5 * ct[0];
            for (size_t j = 1; j < nch; j++) {
                ct[j] = 0.5 * cos(delta * j);
                ct[nc - j] = 0.5 * sin(delta * j);
            }
        }
    }

    static ::std::vector<size_t>& getip(size_t n) {
        size_t j, l, m, pn = 0, nn = n;
        while (nn >>= 1) ++pn;
        if (it.size() <= pn) {
            it.resize(pn+1);
            it.shrink_to_fit();
        }
        ::std::vector<size_t> &ip = it[pn];
        if (ip.size()==0) {
            ip.push_back(0);
            l = n;
            m = 1;
            while ((m << 3) < l) {
                l >>= 1;
                for (j = 0; j < m; j++) {
                    ip.push_back(ip[j] + l);
                }
                m <<= 1;
            }
            ip.push_back(m);
            ip.push_back((m << 3) == l);
            ip.shrink_to_fit();
        }
        return ip;
    }

    static void bitrv2(size_t n, T* a)
    {
        ::std::vector<size_t> const &ip = getip(n);
        size_t j, j1, k, k1, m, m2;
        T xr, xi, yr, yi;
        m = ip[ip.size()-2];
        m2 = 2 * m;
        if (ip[ip.size()-1]) {
            for (k = 0; k < m; k++) {
                for (j = 0; j < k; j++) {
                    j1 = 2 * j + ip[k];
                    k1 = 2 * k + ip[j];
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += m2;
                    k1 += 2 * m2;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += m2;
                    k1 -= m2;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += m2;
                    k1 += 2 * m2;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                }
                j1 = 2 * k + m2 + ip[k];
                k1 = j1 + m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        } else {
            for (k = 1; k < m; k++) {
                for (j = 0; j < k; j++) {
                    j1 = 2 * j + ip[k];
                    k1 = 2 * k + ip[j];
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += m2;
                    k1 += m2;
                    xr = a[j1];
                    xi = a[j1 + 1];
                    yr = a[k1];
                    yi = a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                }
            }
        }
    }

    static void bitrv2(size_t n, const T* a, T* b)
    {
        ::std::vector<size_t> const &ip = getip(n);
        size_t j, j1, k, k1, m, m2;
        m = ip[ip.size()-2];
        m2 = 2 * m;
        if (ip[ip.size()-1]) {
            for (k = 0; k < m; k++) {
                for (j = 0; j < k; j++) {
                    j1 = 2 * j + ip[k];
                    k1 = 2 * k + ip[j];
                    b[k1] = a[j1];
                    b[k1+1] = a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = a[k1 + 1];
                    j1 += m2;
                    k1 += 2 * m2;
                    b[k1] = a[j1];
                    b[k1+1] = a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = a[k1 + 1];
                    j1 += m2;
                    k1 -= m2;
                    b[k1] = a[j1];
                    b[k1+1] = a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = a[k1 + 1];
                    j1 += m2;
                    k1 += 2 * m2;
                    b[k1] = a[j1];
                    b[k1+1] = a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = a[k1 + 1];
                }
                k1 = 2 * k + ip[k];
                b[k1] = a[k1];
                b[k1 + 1] = a[k1 + 1];
                j1 = k1 + m2;
                k1 = j1 + m2;
                b[k1] = a[j1];
                b[k1+1] = a[j1 + 1];
                b[j1] = a[k1];
                b[j1+1] = a[k1 + 1];
                k1 += m2;
                b[k1] = a[k1];
                b[k1 + 1] = a[k1 + 1];
            }
        } else {
            b[0] = a[0];
            b[1] = a[1];
            b[m2] = a[m2];
            b[m2 + 1] = a[m2 + 1];
            for (k = 1; k < m; k++) {
                for (j = 0; j < k; j++) {
                    j1 = 2 * j + ip[k];
                    k1 = 2 * k + ip[j];
                    b[k1] = a[j1];
                    b[k1+1] = a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = a[k1 + 1];
                    j1 += m2;
                    k1 += m2;
                    b[k1] = a[j1];
                    b[k1+1] = a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = a[k1 + 1];
                }
                k1 = 2 * k + ip[k];
                b[k1] = a[k1];
                b[k1 + 1] = a[k1 + 1];
                b[k1 + m2] = a[k1 + m2];
                b[k1 + m2 + 1] = a[k1 + m2 + 1];
            }
        }
    }

    static void bitrv2conj(size_t n, T* a)
    {
        ::std::vector<size_t> const &ip = getip(n);
        size_t j, j1, k, k1, m, m2;
        T xr, xi, yr, yi;
        m = ip[ip.size()-2];
        m2 = 2 * m;
        if (ip[ip.size()-1]) {
            for (k = 0; k < m; k++) {
                for (j = 0; j < k; j++) {
                    j1 = 2 * j + ip[k];
                    k1 = 2 * k + ip[j];
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += m2;
                    k1 += 2 * m2;
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += m2;
                    k1 -= m2;
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += m2;
                    k1 += 2 * m2;
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                }
                k1 = 2 * k + ip[k];
                a[k1 + 1] = -a[k1 + 1];
                j1 = k1 + m2;
                k1 = j1 + m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                k1 += m2;
                a[k1 + 1] = -a[k1 + 1];
            }
        } else {
            a[1] = -a[1];
            a[m2 + 1] = -a[m2 + 1];
            for (k = 1; k < m; k++) {
                for (j = 0; j < k; j++) {
                    j1 = 2 * j + ip[k];
                    k1 = 2 * k + ip[j];
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                    j1 += m2;
                    k1 += m2;
                    xr = a[j1];
                    xi = -a[j1 + 1];
                    yr = a[k1];
                    yi = -a[k1 + 1];
                    a[j1] = yr;
                    a[j1 + 1] = yi;
                    a[k1] = xr;
                    a[k1 + 1] = xi;
                }
                k1 = 2 * k + ip[k];
                a[k1 + 1] = -a[k1 + 1];
                a[k1 + m2 + 1] = -a[k1 + m2 + 1];
            }
        }
    }

    static void bitrv2conj(size_t n, const T* a, T* b)
    {
        ::std::vector<size_t> const &ip = getip(n);
        size_t j, j1, k, k1, m, m2;
        m = ip[ip.size()-2];
        m2 = 2 * m;
        if (ip[ip.size()-1]) {
            for (k = 0; k < m; k++) {
                for (j = 0; j < k; j++) {
                    j1 = 2 * j + ip[k];
                    k1 = 2 * k + ip[j];
                    b[k1] = a[j1];
                    b[k1+1] = -a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = -a[k1 + 1];
                    j1 += m2;
                    k1 += 2 * m2;
                    b[k1] = a[j1];
                    b[k1+1] = -a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = -a[k1 + 1];
                    j1 += m2;
                    k1 -= m2;
                    b[k1] = a[j1];
                    b[k1+1] = -a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = -a[k1 + 1];
                    j1 += m2;
                    k1 += 2 * m2;
                    b[k1] = a[j1];
                    b[k1+1] = -a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = -a[k1 + 1];
                }
                k1 = 2 * k + ip[k];
                b[k1] = a[k1];
                b[k1 + 1] = -a[k1 + 1];
                j1 = k1 + m2;
                k1 = j1 + m2;
                b[k1] = a[j1];
                b[k1+1] = -a[j1 + 1];
                b[j1] = a[k1];
                b[j1+1] = -a[k1 + 1];
                k1 += m2;
                b[k1] = a[k1];
                b[k1 + 1] = -a[k1 + 1];
            }
        } else {
            b[0] = a[0];
            b[1] = -a[1];
            b[m2] = a[m2];
            b[m2 + 1] = -a[m2 + 1];
            for (k = 1; k < m; k++) {
                for (j = 0; j < k; j++) {
                    j1 = 2 * j + ip[k];
                    k1 = 2 * k + ip[j];
                    b[k1] = a[j1];
                    b[k1+1] = -a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = -a[k1 + 1];
                    j1 += m2;
                    k1 += m2;
                    b[k1] = a[j1];
                    b[k1+1] = -a[j1 + 1];
                    b[j1] = a[k1];
                    b[j1+1] = -a[k1 + 1];
                }
                k1 = 2 * k + ip[k];
                b[k1] = a[k1];
                b[k1 + 1] = -a[k1 + 1];
                b[k1 + m2] = a[k1 + m2];
                b[k1 + m2 + 1] = -a[k1 + m2 + 1];
            }
        }
    }

    static void cftfsub(size_t n, T *a)
    {
        size_t j, j1, j2, j3, l;
        T x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

        l = 2;
        if (n > 8) {
            cft1st(n, a);
            l = 8;
            while ((l << 2) < n) {
                cftmdl(n, l, a);
                l <<= 2;
            }
        }
        if ((l << 2) == n) {
            for (j = 0; j < l; j += 2) {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                x0r = a[j] + a[j1];
                x0i = a[j + 1] + a[j1 + 1];
                x1r = a[j] - a[j1];
                x1i = a[j + 1] - a[j1 + 1];
                x2r = a[j2] + a[j3];
                x2i = a[j2 + 1] + a[j3 + 1];
                x3r = a[j2] - a[j3];
                x3i = a[j2 + 1] - a[j3 + 1];
                a[j] = x0r + x2r;
                a[j + 1] = x0i + x2i;
                a[j2] = x0r - x2r;
                a[j2 + 1] = x0i - x2i;
                a[j1] = x1r - x3i;
                a[j1 + 1] = x1i + x3r;
                a[j3] = x1r + x3i;
                a[j3 + 1] = x1i - x3r;
            }
        } else {
            for (j = 0; j < l; j += 2) {
                j1 = j + l;
                x0r = a[j] - a[j1];
                x0i = a[j + 1] - a[j1 + 1];
                a[j] += a[j1];
                a[j + 1] += a[j1 + 1];
                a[j1] = x0r;
                a[j1 + 1] = x0i;
            }
        }
    }

    static void cftbsub(size_t n, T *a)
    {
        size_t j, j1, j2, j3, l;
        T x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

        l = 2;
        if (n > 8) {
            cft1st(n, a);
            l = 8;
            while ((l << 2) < n) {
                cftmdl(n, l, a);
                l <<= 2;
            }
        }
        if ((l << 2) == n) {
            for (j = 0; j < l; j += 2) {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                x0r = a[j] + a[j1];
                x0i = -a[j + 1] - a[j1 + 1];
                x1r = a[j] - a[j1];
                x1i = -a[j + 1] + a[j1 + 1];
                x2r = a[j2] + a[j3];
                x2i = a[j2 + 1] + a[j3 + 1];
                x3r = a[j2] - a[j3];
                x3i = a[j2 + 1] - a[j3 + 1];
                a[j] = x0r + x2r;
                a[j + 1] = x0i - x2i;
                a[j2] = x0r - x2r;
                a[j2 + 1] = x0i + x2i;
                a[j1] = x1r - x3i;
                a[j1 + 1] = x1i - x3r;
                a[j3] = x1r + x3i;
                a[j3 + 1] = x1i + x3r;
            }
        } else {
            for (j = 0; j < l; j += 2) {
                j1 = j + l;
                x0r = a[j] - a[j1];
                x0i = -a[j + 1] + a[j1 + 1];
                a[j] += a[j1];
                a[j + 1] = -a[j + 1] - a[j1 + 1];
                a[j1] = x0r;
                a[j1 + 1] = x0i;
            }
        }
    }

    static void cft1st(size_t n, T *a)
    {
        size_t j, k1, k2;
        T wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
        T x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

        x0r = a[0] + a[2];
        x0i = a[1] + a[3];
        x1r = a[0] - a[2];
        x1i = a[1] - a[3];
        x2r = a[4] + a[6];
        x2i = a[5] + a[7];
        x3r = a[4] - a[6];
        x3i = a[5] - a[7];
        a[0] = x0r + x2r;
        a[1] = x0i + x2i;
        a[4] = x0r - x2r;
        a[5] = x0i - x2i;
        a[2] = x1r - x3i;
        a[3] = x1i + x3r;
        a[6] = x1r + x3i;
        a[7] = x1i - x3r;
        wk1r = wt[2];
        x0r = a[8] + a[10];
        x0i = a[9] + a[11];
        x1r = a[8] - a[10];
        x1i = a[9] - a[11];
        x2r = a[12] + a[14];
        x2i = a[13] + a[15];
        x3r = a[12] - a[14];
        x3i = a[13] - a[15];
        a[8] = x0r + x2r;
        a[9] = x0i + x2i;
        a[12] = x2i - x0i;
        a[13] = x0r - x2r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[10] = wk1r * (x0r - x0i);
        a[11] = wk1r * (x0r + x0i);
        x0r = x3i + x1r;
        x0i = x3r - x1i;
        a[14] = wk1r * (x0i - x0r);
        a[15] = wk1r * (x0i + x0r);
        k1 = 0;
        for (j = 16; j < n; j += 16) {
            k1 += 2;
            k2 = 2 * k1;
            wk2r = wt[k1];
            wk2i = wt[k1 + 1];
            wk1r = wt[k2];
            wk1i = wt[k2 + 1];
            wk3r = wk1r - 2 * wk2i * wk1i;
            wk3i = 2 * wk2i * wk1r - wk1i;
            x0r = a[j] + a[j + 2];
            x0i = a[j + 1] + a[j + 3];
            x1r = a[j] - a[j + 2];
            x1i = a[j + 1] - a[j + 3];
            x2r = a[j + 4] + a[j + 6];
            x2i = a[j + 5] + a[j + 7];
            x3r = a[j + 4] - a[j + 6];
            x3i = a[j + 5] - a[j + 7];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j + 4] = wk2r * x0r - wk2i * x0i;
            a[j + 5] = wk2r * x0i + wk2i * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j + 2] = wk1r * x0r - wk1i * x0i;
            a[j + 3] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j + 6] = wk3r * x0r - wk3i * x0i;
            a[j + 7] = wk3r * x0i + wk3i * x0r;
            wk1r = wt[k2 + 2];
            wk1i = wt[k2 + 3];
            wk3r = wk1r - 2 * wk2r * wk1i;
            wk3i = 2 * wk2r * wk1r - wk1i;
            x0r = a[j + 8] + a[j + 10];
            x0i = a[j + 9] + a[j + 11];
            x1r = a[j + 8] - a[j + 10];
            x1i = a[j + 9] - a[j + 11];
            x2r = a[j + 12] + a[j + 14];
            x2i = a[j + 13] + a[j + 15];
            x3r = a[j + 12] - a[j + 14];
            x3i = a[j + 13] - a[j + 15];
            a[j + 8] = x0r + x2r;
            a[j + 9] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j + 12] = -wk2i * x0r - wk2r * x0i;
            a[j + 13] = -wk2i * x0i + wk2r * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j + 10] = wk1r * x0r - wk1i * x0i;
            a[j + 11] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j + 14] = wk3r * x0r - wk3i * x0i;
            a[j + 15] = wk3r * x0i + wk3i * x0r;
        }
    }

    static void cftmdl(size_t n, size_t l, T *a)
    {
        size_t j, j1, j2, j3, k, k1, k2, m, m2;
        T wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
        T x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

        m = l << 2;
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i - x2i;
            a[j1] = x1r - x3i;
            a[j1 + 1] = x1i + x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
        }
        wk1r = wt[2];
        for (j = m; j < l + m; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            a[j2] = x2i - x0i;
            a[j2 + 1] = x0r - x2r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * (x0r - x0i);
            a[j1 + 1] = wk1r * (x0r + x0i);
            x0r = x3i + x1r;
            x0i = x3r - x1i;
            a[j3] = wk1r * (x0i - x0r);
            a[j3 + 1] = wk1r * (x0i + x0r);
        }
        k1 = 0;
        m2 = 2 * m;
        for (k = m2; k < n; k += m2) {
            k1 += 2;
            k2 = 2 * k1;
            wk2r = wt[k1];
            wk2i = wt[k1 + 1];
            wk1r = wt[k2];
            wk1i = wt[k2 + 1];
            wk3r = wk1r - 2 * wk2i * wk1i;
            wk3i = 2 * wk2i * wk1r - wk1i;
            for (j = k; j < l + k; j += 2) {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                x0r = a[j] + a[j1];
                x0i = a[j + 1] + a[j1 + 1];
                x1r = a[j] - a[j1];
                x1i = a[j + 1] - a[j1 + 1];
                x2r = a[j2] + a[j3];
                x2i = a[j2 + 1] + a[j3 + 1];
                x3r = a[j2] - a[j3];
                x3i = a[j2 + 1] - a[j3 + 1];
                a[j] = x0r + x2r;
                a[j + 1] = x0i + x2i;
                x0r -= x2r;
                x0i -= x2i;
                a[j2] = wk2r * x0r - wk2i * x0i;
                a[j2 + 1] = wk2r * x0i + wk2i * x0r;
                x0r = x1r - x3i;
                x0i = x1i + x3r;
                a[j1] = wk1r * x0r - wk1i * x0i;
                a[j1 + 1] = wk1r * x0i + wk1i * x0r;
                x0r = x1r + x3i;
                x0i = x1i - x3r;
                a[j3] = wk3r * x0r - wk3i * x0i;
                a[j3 + 1] = wk3r * x0i + wk3i * x0r;
            }
            wk1r = wt[k2 + 2];
            wk1i = wt[k2 + 3];
            wk3r = wk1r - 2 * wk2r * wk1i;
            wk3i = 2 * wk2r * wk1r - wk1i;
            for (j = k + m; j < l + (k + m); j += 2) {
                j1 = j + l;
                j2 = j1 + l;
                j3 = j2 + l;
                x0r = a[j] + a[j1];
                x0i = a[j + 1] + a[j1 + 1];
                x1r = a[j] - a[j1];
                x1i = a[j + 1] - a[j1 + 1];
                x2r = a[j2] + a[j3];
                x2i = a[j2 + 1] + a[j3 + 1];
                x3r = a[j2] - a[j3];
                x3i = a[j2 + 1] - a[j3 + 1];
                a[j] = x0r + x2r;
                a[j + 1] = x0i + x2i;
                x0r -= x2r;
                x0i -= x2i;
                a[j2] = -wk2i * x0r - wk2r * x0i;
                a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
                x0r = x1r - x3i;
                x0i = x1i + x3r;
                a[j1] = wk1r * x0r - wk1i * x0i;
                a[j1 + 1] = wk1r * x0i + wk1i * x0r;
                x0r = x1r + x3i;
                x0i = x1i - x3r;
                a[j3] = wk3r * x0r - wk3i * x0i;
                a[j3 + 1] = wk3r * x0i + wk3i * x0r;
            }
        }
    }

    static void rftfsub(size_t n, T *a, size_t nc)
    {
        size_t j, k, kk, ks, m;
        T wkr, wki, xr, xi, yr, yi;

        m = n >> 1;
        ks = 2 * nc / m;
        kk = 0;
        for (j = 2; j < m; j += 2) {
            k = n - j;
            kk += ks;
            wkr = 0.5 - ct[nc - kk];
            wki = ct[kk];
            xr = a[j] - a[k];
            xi = a[j + 1] + a[k + 1];
            yr = wkr * xr - wki * xi;
            yi = wkr * xi + wki * xr;
            a[j] -= yr;
            a[j + 1] -= yi;
            a[k] += yr;
            a[k + 1] -= yi;
        }
    }

    static void rftbsub(size_t n, T *a, size_t nc)
    {
        size_t j, k, kk, ks, m;
        T wkr, wki, xr, xi, yr, yi;

        a[1] = -a[1];
        m = n >> 1;
        ks = 2 * nc / m;
        kk = 0;
        for (j = 2; j < m; j += 2) {
            k = n - j;
            kk += ks;
            wkr = 0.5 - ct[nc - kk];
            wki = ct[kk];
            xr = a[j] - a[k];
            xi = a[j + 1] + a[k + 1];
            yr = wkr * xr + wki * xi;
            yi = wkr * xi - wki * xr;
            a[j] -= yr;
            a[j + 1] = yi - a[j + 1];
            a[k] += yr;
            a[k + 1] = yi - a[k + 1];
        }
        a[m + 1] = -a[m + 1];
    }

    static void rftbsub(size_t n, const T *a, T *b, size_t nc)
    {
        size_t j, k, kk, ks, m;
        T wkr, wki, xr, xi, yr, yi;

        // [0] and [1] are copied by caller
        b[1] = -b[1];
        m = n >> 1;
        ks = 2 * nc / m;
        kk = 0;
        for (j = 2; j < m; j += 2) {
            k = n - j;
            kk += ks;
            wkr = 0.5 - ct[nc - kk];
            wki = ct[kk];
            xr = a[j] - a[k];
            xi = a[j + 1] + a[k + 1];
            yr = wkr * xr + wki * xi;
            yi = wkr * xi - wki * xr;
            b[j] = a[j] - yr;
            b[j + 1] = yi - a[j + 1];
            b[k] = a[k] + yr;
            b[k + 1] = yi - a[k + 1];
        }
        b[m] = a[m];
        b[m + 1] = -a[m + 1];
    }

    static void dctsub(size_t n, T *a, size_t nc)
    {
        size_t j, k, kk, ks, m;
        T wkr, wki, xr;

        m = n >> 1;
        ks = nc / n;
        kk = 0;
        for (j = 1; j < m; j++) {
            k = n - j;
            kk += ks;
            wkr = ct[kk] - ct[nc - kk];
            wki = ct[kk] + ct[nc - kk];
            xr = wki * a[j] - wkr * a[k];
            a[j] = wkr * a[j] + wki * a[k];
            a[k] = xr;
        }
        a[m] *= ct[0];
    }

    static void dstsub(size_t n, T *a, size_t nc)
    {
        size_t j, k, kk, ks, m;
        T wkr, wki, xr;

        m = n >> 1;
        ks = nc / n;
        kk = 0;
        for (j = 1; j < m; j++) {
            k = n - j;
            kk += ks;
            wkr = ct[kk] - ct[nc - kk];
            wki = ct[kk] + ct[nc - kk];
            xr = wki * a[k] - wkr * a[j];
            a[k] = wkr * a[k] + wki * a[j];
            a[j] = xr;
        }
        a[m] *= ct[0];
    }
};

static ::std::vector<::std::vector<size_t>> FFT_Impl_it;
template <typename T>
::std::vector<::std::vector<size_t>> &FFT_Impl<T>::it = FFT_Impl_it;
template <typename T>
size_t FFT_Impl<T>::wt_n = 0;
template <typename T>
::std::vector<T> FFT_Impl<T>::wt;
template <typename T>
size_t FFT_Impl<T>::ct_n = 0;
template <typename T>
::std::vector<T> FFT_Impl<T>::ct;

/// @endcond DSPP_IMPL

/// @brief Fast Fourier/Cosine/Sine Transforms
/// @details
/// Performs a transform on arrays with power-of-two sizes.
/// Uses Sande-Tukey decimation in frequency algorithm.
///
/// Example usage:
///
///     std::vector<std::complex<float>> data(512);
///     dspp::FFT::dft(data);
///
/// C-style arrays are supported:
///
///     std::complex<float> data[512];
///     dspp::FFT::dft(data, 1); // 1 = inverse DFT
///
/// Out-of-place transform:
///
///     std::array<double,64> in;
///     std::array<std::complex<double>,32> out;
///     dspp::FFT::dft(in, out); // real transform
///
/// Planning is automatically performed as necessary. The first transform
/// of a specific size and type will take longer than subsequent transforms.
///
/// Each of the transform types has an unsafe version. If possible, always use
/// the convenience overloads which engage the type safety of your compiler.
namespace FFT {

/// Complex discrete Fourier transform. In-place unsafe version.
template<typename T>
inline void dft(size_t n, ::std::complex<T> *data, int isgn = -1)
{
    FFT_Impl<T>::cdft(n*2, isgn, (T*)data);
}

/// Complex discrete Fourier transform. Out-of-place unsafe version.
template<typename T>
inline void dft(size_t insize, const ::std::complex<T> *in, ::std::complex<T> *out, int isgn = -1)
{
    FFT_Impl<T>::cdft(insize*2, isgn, (T*)in, (T*)out);
}

/// @brief Real discrete Fourier transform. In-place unsafe version.
/// @details Value for data[n].real() is packed in data[0].imag().
/// @param isgn -1 is real to complex.<br>
///             1 is complex to real.
template<typename T>
inline void dft(size_t n, T *data, int isgn = -1)
{
    FFT_Impl<T>::rdft(n, -isgn, data);
}

/// Real to complex discrete Fourier transform. Out-of-place unsafe version.
template<typename T>
inline void dft(size_t insize, const T *in, ::std::complex<T> *out, int unused = 0)
{
    FFT_Impl<T>::rdft(insize, 1, in, (T*)out);
}

/// Complex to real discrete Fourier transform. Out-of-place unsafe version.
template<typename T>
inline void dft(size_t insize, const ::std::complex<T> *in, T *out, int unused = 0)
{
    FFT_Impl<T>::rdft(insize*2, -1, (T*)in, out);
}

/// Discrete Fourier transform. In-place for C-style array.
template<typename T, size_t N>
inline void dft(T (&data)[N], int isgn = -1) {
    dft(N, (T*)&data, isgn);
}

/// Discrete Fourier transform. In-place for STL container.
template<typename T>
inline void dft(T &data, int isgn = -1)
{
    dft(data.size(), data.data(), isgn);
}

/// Discrete Fourier transform. Out-of-place for C-style arrays.
template<typename T1, size_t N1, typename T2, size_t N2>
inline void dft(const T1 (&in)[N1], T2 (&out)[N2], int isgn = -1) {
    assert(sizeof(in) == sizeof(out));
    dft(N1, (T1*)&in, (T2*)&out, isgn);
}

/// Discrete Fourier transform. Out-of-place for STL containers.
template<typename T1, typename T2>
inline void dft(const T1 &in, T2 &out, int isgn = -1)
{
    assert(in.size()*sizeof(*in.data()) == out.size()*sizeof(*out.data()));
    dft(in.size(), in.data(), out.data(), isgn);
}

/// Discrete cosine transform. In-place unsafe version.
template<typename T>
inline void dct(size_t n, T *data, int isgn = -1)
{
    FFT_Impl<T>::ddct(n, isgn, data);
}

/// Discrete sine transform. In-place unsafe version.
template<typename T>
inline void dst(size_t n, T *data, int isgn = -1)
{
    FFT_Impl<T>::ddst(n, isgn, data);
}


} /* namespace FFT */
} /* namespace dspp */

#endif /* DSPP_FFT_HPP */
