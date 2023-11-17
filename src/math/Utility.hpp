#pragma once

#include "rack.hpp"

using rack::simd::float_4;

// https://web.archive.org/web/20200628195036/http://mooooo.ooo/chebyshev-sine-approximation/
// Thanks to Colin Wallace for doing some great math
template<typename T>
T sin2pi_chebyshev(T x) {
    const T CHEBYSHEV_COEFFS[6] = {
        T(-3.1415926444234477f),   // x
        T(2.0261194642649887f),    // x^3
        T(-0.5240361513980939f),   // x^5
        T(0.0751872634325299f),    // x^7
        T(-0.006860187425683514f), // x^9
        T(0.000385937753182769f),  // x^11
    };

    x = (-x * 2.f) + 1.f;
    auto x2 = x * x;

    // Perfectly fine without the x^9 or x^11 terms
    // Error is at most 0.000482 at values 0.0482 from the boundaries 0 and 1
    auto p7 = CHEBYSHEV_COEFFS[3];
    auto p5 = p7 * x2 + CHEBYSHEV_COEFFS[2];
    auto p3 = p5 * x2 + CHEBYSHEV_COEFFS[1];
    auto p1 = p3 * x2 + CHEBYSHEV_COEFFS[0];
    return (x - 1.f) * (x + 1.f) * p1 * x;
}

template<typename T>
void sincos2pi_chebyshev(T x, T &sinOut, T &cosOut) {
    const T SIN_COEFFS[4] = {
        T(-3.1415926444234477f),  // x
        T(2.0261194642649887f),   // x^3
        T(-0.5240361513980939f),  // x^5
        T(0.0751872634325299f),   // x^7
    };

    const T COS_COEFFS[4] = {
        T(-0.318309887112536f),   // overall scale = 1/SIN_COEFFS[0]
        T(4.052238928529978f),    // x^2
        T(-2.096144605592376f),   // x^4
        T(0.451123580595179f),    // x^6
    };

    x = (-x * 2.f) + 1.f;
    auto x2 = x * x;

    // Sin
    // Error is at most 0.000482 at values 0.0482 from the boundaries 0 and 1
    auto p7 = SIN_COEFFS[3];
    auto p5 = p7 * x2 + SIN_COEFFS[2];
    auto p3 = p5 * x2 + SIN_COEFFS[1];
    auto p1 = p3 * x2 + SIN_COEFFS[0];
    auto sinMult = x2 * x - x; // (x - 1) * (x + 1) * x
    sinOut = p1 * sinMult;

    // Cos (product rule)
    // Error is at most 0.004122 at the boundaries 0 and 1
    auto p6 = COS_COEFFS[3];
    auto p4 = p6 * x2 + COS_COEFFS[2];
    auto p2 = p4 * x2 + COS_COEFFS[1];
    auto deriv = (sinMult * x * p2) + (3.f * x2 - 1.f) * p1;
    cosOut = COS_COEFFS[0] * deriv;
}

// See: https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction
float sum_float4(const float_4 &x) {
#ifdef __SSE3__
    __m128 shuf = _mm_movehdup_ps(x.v);
#else
    __m128 shuf = _mm_shuffle_ps(x.v, x.v, _MM_SHUFFLE(2, 3, 0, 1));
#endif
    __m128 sums = _mm_add_ps(x.v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums);
    sums        = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}
