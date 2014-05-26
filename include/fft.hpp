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

namespace dspp {
    namespace FFT {

        /// Computes and caches twiddle factors.
        template<typename T, size_t N, int D>
        class Twiddle4 {
        public:
            static std::array<std::complex<T>, N> t1;
            static std::array<std::complex<T>, N> t2;
            static std::array<std::complex<T>, N> t3;
            Twiddle4() {
                // Initialize values on first instance.
                if (!t1[0].real()) {
                    T theta = M_PI*2*D/(N*4);
                    for (size_t i=0; i < N; ++i) {
                        T phi = i*theta;
                        t1[i] = std::complex<T>(cos(phi), sin(phi));
                        t2[i] = std::complex<T>(cos(phi*2), sin(phi*2));
                        t3[i] = std::complex<T>(cos(phi*3), sin(phi*3));
                    }
                }
            }
        };

        // Storage for twiddle factors
        template<typename T, size_t N, int D>
        std::array<std::complex<T>, N> Twiddle4<T, N, D>::t1 {{0}};
        template<typename T, size_t N, int D>
        std::array<std::complex<T>, N> Twiddle4<T, N, D>::t2;
        template<typename T, size_t N, int D>
        std::array<std::complex<T>, N> Twiddle4<T, N, D>::t3;

        /// Recursive template for mixing.
        template<typename T, size_t N, int D>
        class Radix4 {
            static const size_t N4 = N/4;
            const Twiddle4<T, N4, D> t;
            Radix4<T, N4, D> next;
            /// Limited range (fast) multiplication of complex numbers.
            static inline std::complex<T> multiply(const std::complex<T>& z, const std::complex<T>& w)
            {
                T a = z.real();
                T b = z.imag();
                T c = w.real();
                T d = w.imag();
                return std::complex<T>(a*c - b*d, a*d + b*c);
            }
            /// Simplified multiplication for direction product.
            static inline std::complex<T> direction(const std::complex<T>& z)
            {
                if (D>0) return std::complex<T>(-z.imag(), z.real());
                else return std::complex<T>(z.imag(), -z.real());
            }
        public:
            inline void mix(std::complex<T>* data) {
                size_t i1 = N4;
                size_t i2 = N4 + N4;
                size_t i3 = i2 + N4;
                next.mix(data);
                next.mix(data+i1);
                next.mix(data+i2);
                next.mix(data+i3);
                // Index 0 twiddles are always (1+0i).
                std::complex<T> a0 = data[0];
                std::complex<T> a1 = data[i2];
                std::complex<T> a2 = data[i1];
                std::complex<T> a3 = data[i3];
                std::complex<T> b0 = (a1+a3);
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
                    a1 = multiply(data[i2], t.t1[i0]);
                    a2 = multiply(data[i1], t.t2[i0]);
                    a3 = multiply(data[i3], t.t3[i0]);
                    b0 = (a1+a3);
                    b1 = direction(a1-a3);
                    data[i0] = a0 + a2 + b0;
                    data[i1] = a0 - a2 + b1;
                    data[i2] = a0 + a2 - b0;
                    data[i3] = a0 - a2 - b1;
                }
            }
        };

        /// Terminates template recursion for powers of 2.
        template<typename T, int D>
        class Radix4<T, 2, D> {
        public:
            inline void mix(std::complex<T>* data) {
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
            void mix(std::complex<T>* data) {
                // Do nothing.
            }
        };

        /// Fast Fourier Transform.
        template<typename T, size_t N>
        class FFT {
            static_assert((N > 1) & !(N & (N - 1)), "Array size must be a power of two.");
            static void reindex(std::array<std::complex<T>, N> &data) {
                for (size_t i=0, j=0; i<N; ++i) {
                    if (j>i) {
                        std::swap(data[i], data[j]);
                    }
                    size_t m = N>>1;
                    while (m>=1 && j>=m) {
                        j -= m;
                        m >>= 1;
                    }
                    j += m;
                }
            }
            static void reindex(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) {
                for (size_t i=0, j=0; i<N; ++i) {
                    out[i] = in[j];
                    size_t m = N>>1;
                    while (m>=1 && j>=m) {
                        j -= m;
                        m >>= 1;
                    }
                    j += m;
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
            virtual void plan_dft() {
                Radix4<T, N, -1>();
            }
            virtual void plan_idft() {
                Radix4<T, N, 1>();
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
                Radix4<T, N, -1> mixer;
                reindex(data);
                mixer.mix(&data[0]);
            }
            static void idft(std::array<std::complex<T>, N> &data) {
                Radix4<T, N, 1> mixer;
                reindex(data);
                mixer.mix(&data[0]);
            }
            static void dft(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) {
                Radix4<T, N, -1> mixer;
                reindex(in, out);
                mixer.mix(&out[0]);
            }
            static void idft(const std::array<std::complex<T>, N> &in, std::array<std::complex<T>, N> &out) {
                Radix4<T, N, 1> mixer;
                reindex(in, out);
                mixer.mix(&out[0]);
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
