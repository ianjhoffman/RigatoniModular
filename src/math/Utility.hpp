#pragma once

#include "rack.hpp"

using rack::simd::float_4;

// https://web.archive.org/web/20200628195036/http://mooooo.ooo/chebyshev-sine-approximation/
// Thanks to Colin Wallace for doing some great math
const float_4 CHEBYSHEV_COEFFS[6] = {
	float_4(-3.1415926444234477f),   // x
	float_4(2.0261194642649887f),    // x^3
	float_4(-0.5240361513980939f),   // x^5
	float_4(0.0751872634325299f),    // x^7
	float_4(-0.006860187425683514f), // x^9
	float_4(0.000385937753182769f),  // x^11
};

float_4 sin2pi_chebyshev(float_4 x) {
	x = (-x * 2.f) + 1.f;
	auto x2 = x * x;
    auto p11 = CHEBYSHEV_COEFFS[5];
    auto p9 = p11 * x2 + CHEBYSHEV_COEFFS[4];
    auto p7 = p9 * x2  + CHEBYSHEV_COEFFS[3];
    auto p5 = p7 * x2  + CHEBYSHEV_COEFFS[2];
    auto p3 = p5 * x2  + CHEBYSHEV_COEFFS[1];
    auto p1 = p3 * x2  + CHEBYSHEV_COEFFS[0];
    return (x - 1.f) * (x + 1.f) * p1 * x;
}

inline float sum_float4(float_4 x) {
	return x[0] + x[1] + x[2] + x[3];
}