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

#ifndef dspp_fft_hpp
#define dspp_fft_hpp

#include <complex>
#include <array>
#include <cmath>
#include "util.hpp"

namespace dspp {
    namespace FFT {

        /// Recursive template for butterfly mixing.
        template<typename T, size_t N, int D>
        class Radix4 {
            static const std::array<std::complex<T>, N/4> t1;
            static const std::array<std::complex<T>, N/4> t2;
            static const std::array<std::complex<T>, N/4> t3;
            static const size_t N4 = N/4;
            Radix4<T, N4, D> next;
            /// Simplified multiplication for direction product.
            static inline std::complex<T> direction(const std::complex<T>& z)
            {
                if (D>0) return std::complex<T>(-z.imag(), z.real());
                else return std::complex<T>(z.imag(), -z.real());
            }
        public:
            void operator()(std::complex<T>* data) {
                size_t i1 = N4;
                size_t i2 = N4 + N4;
                size_t i3 = i2 + N4;
                next(data);
                next(data+i1);
                next(data+i2);
                next(data+i3);
                // Index 0 twiddles are always (1+0i).
                std::complex<T> a0 = data[0];
                std::complex<T> a2 = data[i1];
                std::complex<T> a1 = data[i2];
                std::complex<T> a3 = data[i3];
                std::complex<T> b0 = a1 + a3;
                std::complex<T> b1 = direction(a1-a3);
                data[0] = a0 + a2 + b0;
                data[i1] = a0 - a2 + b1;
                data[i2] = a0 + a2 - b0;
                data[i3] = a0 - a2 - b1;
                // Index 1+ must multiply twiddles.
                for (size_t i0=1; i0 < N4; ++i0) {
                    i1 = i0 + N4;
                    i2 = i1 + N4;
                    i3 = i2 + N4;
                    a0 = data[i0];
                    a2 = mul(data[i1], t2[i0]);
                    a1 = mul(data[i2], t1[i0]);
                    a3 = mul(data[i3], t3[i0]);
                    b0 = a1 + a3;
                    b1 = direction(a1-a3);
                    data[i0] = a0 + a2 + b0;
                    data[i1] = a0 - a2 + b1;
                    data[i2] = a0 + a2 - b0;
                    data[i3] = a0 - a2 - b1;
                }
            }
        };

        /// Computes twiddle factors
        template<typename T>
        static std::vector<std::complex<T>> twiddle4(double a, size_t n, int d) {
            std::vector<std::complex<T>> twids(n/4);
            double theta = two_pi<double>()*d/n;
            for (size_t i=0; i < n/4; ++i) {
                double phi = theta * a * i;
                twids[i] = std::complex<T>(cos(phi), sin(phi));
            }
            return twids;
        }

        // Storage and initialization for twiddle factors
        template<typename T, size_t N, int D>
        const std::array<std::complex<T>, N/4> Radix4<T, N, D>::t1(
            *reinterpret_cast<std::array<std::complex<T>, N/4>*>(twiddle4<T>(1,N,D).data())
        );
        template<typename T, size_t N, int D>
        const std::array<std::complex<T>, N/4> Radix4<T, N, D>::t2(
            *reinterpret_cast<std::array<std::complex<T>, N/4>*>(twiddle4<T>(2,N,D).data())
        );
        template<typename T, size_t N, int D>
        const std::array<std::complex<T>, N/4> Radix4<T, N, D>::t3(
            *reinterpret_cast<std::array<std::complex<T>, N/4>*>(twiddle4<T>(3,N,D).data())
        );

        /// Terminates template recursion when not power of 4.
        template<typename T, int D>
        class Radix4<T, 2, D> {
        public:
            inline void operator()(std::complex<T>* data) {
                std::complex<T> a0 = data[0];
                std::complex<T> a1 = data[1];
                data[0] = a0 + a1;
                data[1] = a0 - a1;
            }
        };

        /// Terminates template recursion for powers of 4.
        template<typename T, int D>
        class Radix4<T, 1, D> {
        public:
            inline void operator()(std::complex<T>* data) {
                // Do nothing.
            }
        };

        // Computes the bit reversal pattern for size N
        template<size_t N>
        struct BitReverse {
            static std::vector<size_t> pattern;
            static bool ispow4;
            static std::vector<size_t> calc() {
                std::vector<size_t> l;
                l.push_back(0);
                size_t m = 1;
                size_t n = N;
                while ((m << 2) < n) {
                    n >>= 1;
                    for (size_t j = 0; j < m; j++) {
                        l.push_back(l[j] + n);
                    }
                    m <<= 1;
                }
                ispow4 = (m << 2) == n;
                l.shrink_to_fit();
                return l;
            }
        };

        template<size_t N>
        std::vector<size_t> BitReverse<N>::pattern = BitReverse<N>::calc();
        template<size_t N>
        bool BitReverse<N>::ispow4;

        /// Fast Fourier Transform.
        template<typename T, size_t N>
        class FFT : public BitReverse<N> {
            static_assert((N > 1) & !(N & (N - 1)), "Array size must be a power of two.");
            using BitReverse<N>::pattern;
            using BitReverse<N>::ispow4;

            static void reindex(std::array<std::complex<T>, N> &data) {
                size_t j1, k1, m = pattern.size();
                if (ispow4) {
                    size_t m2 = 2 * m;
                    for (size_t k = 0; k < m; k++) {
                        for (size_t j = 0; j < k; j++) {
                            j1 = j + pattern[k];
                            k1 = k + pattern[j];
                            std::swap(data[j1], data[k1]);
                            j1 += m;
                            k1 += m2;
                            std::swap(data[j1], data[k1]);
                            j1 += m;
                            k1 -= m;
                            std::swap(data[j1], data[k1]);
                            j1 += m;
                            k1 += m2;
                            std::swap(data[j1], data[k1]);
                        }
                        j1 = k + m + pattern[k];
                        k1 = j1 + m;
                        std::swap(data[j1], data[k1]);
                    }
                } else {
                    for (size_t k = 1; k < m; k++) {
                        for (size_t j = 0; j < k; j++) {
                            j1 = j + pattern[k];
                            k1 = k + pattern[j];
                            std::swap(data[j1], data[k1]);
                            j1 += m;
                            k1 += m;
                            std::swap(data[j1], data[k1]);
                        }
                    }
                }
            }

            static void reindex(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) {
                size_t j1, k1, m = pattern.size();
                if (ispow4) {
                    size_t m2 = 2 * m;
                    for (size_t k = 0; k < m; k++) {
                        for (size_t j = 0; j < k; j++) {
                            j1 = j + pattern[k];
                            k1 = k + pattern[j];
                            out[j1] = in[k1];
                            out[k1] = in[j1];
                            j1 += m;
                            k1 += m2;
                            out[j1] = in[k1];
                            out[k1] = in[j1];
                            j1 += m;
                            k1 -= m;
                            out[j1] = in[k1];
                            out[k1] = in[j1];
                            j1 += m;
                            k1 += m2;
                            out[j1] = in[k1];
                            out[k1] = in[j1];
                        }
                        k1 = k + pattern[k];
                        out[k1] = in[k1];
                        j1 = k1 + m;
                        k1 = j1 + m;
                        out[j1] = in[k1];
                        out[k1] = in[j1];
                        k1 += m;
                        out[k1] = in[k1];
                    }
                } else {
                    out[0] = in[0];
                    out[m] = in[m];
                    for (size_t k = 1; k < m; k++) {
                        for (size_t j = 0; j < k; j++) {
                            j1 = j + pattern[k];
                            k1 = k + pattern[j];
                            out[j1] = in[k1];
                            out[k1] = in[j1];
                            j1 += m;
                            k1 += m;
                            out[j1] = in[k1];
                            out[k1] = in[j1];
                        }
                        k1 = k + pattern[k];
                        out[k1] = in[k1];
                        out[k1+m] = in[k1+m];
                    }
                }
            }

            const std::array<std::complex<T>, N> *in;
            std::array<std::complex<T>, N> *out;
        public:
            FFT(std::array<std::complex<T>, N> &data) :
                in(nullptr), out(&data) {
            }
            FFT(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) :
                in(&in), out(&out) {
            }
            virtual void dft() {
                if (in==nullptr) dft(*out);
                else dft(*in, *out);
            }
            virtual void idft() {
                if (in==nullptr) idft(*out);
                else idft(*in, *out);
            }
            static void dft(std::array<std::complex<T>, N> &data) {
                Radix4<T, N, -1> butterfly;
                reindex(data);
                butterfly(&data[0]);
            }
            static void idft(std::array<std::complex<T>, N> &data) {
                Radix4<T, N, 1> butterfly;
                reindex(data);
                butterfly(&data[0]);
            }
            static void dft(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) {
                Radix4<T, N, -1> butterfly;
                reindex(in, out);
                butterfly(&out[0]);
            }
            static void idft(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) {
                Radix4<T, N, 1> butterfly;
                reindex(in, out);
                butterfly(&out[0]);
            }
        };

        /// Discrete Fourier transform.
        template<typename T, size_t N>
        inline void dft(std::array<std::complex<T>, N> &data) {
            FFT<T, N>::dft(data);
        }

        /// Discrete Fourier transform.
        template<typename T, size_t N>
        inline void dft(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) {
            FFT<T, N>::dft(in, out);
        }

        /// Discrete Fourier transform.
        template<typename T, size_t N>
        inline void dft(std::complex<T> (&data)[N]) {
            FFT<T, N>::dft(*reinterpret_cast<std::array<std::complex<T>, N>*>(&data));
        }

        /// Discrete Fourier transform.
        template<typename T, size_t N>
        inline void dft(const std::complex<T> (&in)[N], std::complex<T> (&out)[N]) {
            FFT<T, N>::dft(*reinterpret_cast<std::array<std::complex<T>, N>*>(&in), *reinterpret_cast<std::array<std::complex<T>, N>*>(&out));
        }

        /// Inverse discrete Fourier transform.
        template<typename T, size_t N>
        inline void idft(std::array<std::complex<T>, N> &data) {
            FFT<T, N>::idft(data);
        }

        /// Inverse discrete Fourier transform.
        template<typename T, size_t N>
        inline void idft(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) {
            FFT<T, N>::idft(in, out);
        }

        /// Inverse discrete Fourier transform.
        template<typename T, size_t N>
        inline void idft(std::complex<T> (&data)[N]) {
            FFT<T, N>::idft(*reinterpret_cast<std::array<std::complex<T>, N>*>(&data));
        }

        /// Inverse discrete Fourier transform.
        template<typename T, size_t N>
        inline void idft(const std::complex<T> (&in)[N], std::complex<T> (&out)[N]) {
            FFT<T, N>::idft(*reinterpret_cast<std::array<std::complex<T>, N>*>(&in), *reinterpret_cast<std::array<std::complex<T>, N>*>(&out));
        }

    }
}
#endif
