//
//  Copyright (c) 2014 David Turnbull AE9RB. All rights reserved.
//

#ifndef DSPP_MIXER_HPP
#define DSPP_MIXER_HPP

#include <cmath>
#include "util.hpp"

namespace dspp {
    
template<class Container>
struct Mixer {
    typedef typename Container::value_type cx_type;
    typedef typename cx_type::value_type fp_type;
private:
    // FIXUP_RATE is how often we correct the NCO for rounding errors
    static const int FIXUP_RATE = (1 << 4) - 1;
    fp_type _rate;
    fp_type _freq;
    cx_type nco;
    cx_type clk;
    
    void
    compute_clk() {
        fp_type inc = two_pi<fp_type>() * _freq / _rate;
        clk = cx_type(cos(inc), sin(inc));
    }
    
public:
    Mixer(fp_type rate = 96000, fp_type freq = 0) :
    _rate(rate),
    _freq(freq),
    nco(cx_type(1.0, 0.0))
    {
        compute_clk();
    }
    
    void
    operator()(Container& data) {
        auto fixup_counter = FIXUP_RATE;
        for (auto& d : data) {
            if (!(++fixup_counter & FIXUP_RATE)) {
                fp_type gain = 2.0 - (nco.real()*nco.real() + nco.imag()*nco.imag());
                nco = cx_type(nco.real()*gain, nco.imag()*gain);
            }
            nco = mul(nco, clk);
            d = mul(d, nco);
        }
    }
    
    fp_type
    rate(const fp_type rate) {
        _rate = rate;
        compute_clk();
        return _rate;
    }
    
    fp_type
    rate() const {
        return _rate;
    }
    
    fp_type
    freq(const fp_type freq) {
        _freq = freq;
        compute_clk();
        return _freq;
    }
    
    fp_type
    freq() const {
        return _freq;
    }
};

}

#endif


