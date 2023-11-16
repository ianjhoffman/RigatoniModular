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

    // auto p11 = CHEBYSHEV_COEFFS[5];
    // auto p9 = p11 * x2 + CHEBYSHEV_COEFFS[4];
	// auto p7 = p9 * x2  + CHEBYSHEV_COEFFS[3];
	auto p7 = CHEBYSHEV_COEFFS[3];
    auto p5 = p7 * x2  + CHEBYSHEV_COEFFS[2];
    auto p3 = p5 * x2  + CHEBYSHEV_COEFFS[1];
    auto p1 = p3 * x2  + CHEBYSHEV_COEFFS[0];
    return (x - 1.f) * (x + 1.f) * p1 * x;
}

template<typename T>
T sinToCos(T sinPhase, T sinVal) {
	auto toSqrt = T(1) - (sinVal * sinVal);
	auto cosRectified = toSqrt * simd::rsqrt(toSqrt);
	return simd::ifelse((sinPhase - T(0.25f)) <= T(0.5f), -cosRectified, cosRectified);
}

inline float sum_float4(float_4 x) {
	return x[0] + x[1] + x[2] + x[3];
}

template<int SIZE>
struct SinTable {
	static constexpr float FRAC = 1.f / SIZE;
	static constexpr int TABLE_SIZE = SIZE + 1;
	static constexpr int QUAD = SIZE / 4;

	float table[TABLE_SIZE];

	SinTable() {
		for (int i = 0; i < TABLE_SIZE; i++) {
			table[i] = sin2pi_chebyshev(i * FRAC);
		}
	}

	inline void sinCosBounds(
		int32_t integral, float *leftSin, float *rightSin, float *leftCos, float *rightCos
	) {
		*leftSin = table[integral];
		*rightSin = table[(integral + 1) % TABLE_SIZE];
		*leftCos = table[(integral + QUAD + 1) % TABLE_SIZE];
		*rightCos = table[(integral + QUAD + 1) % TABLE_SIZE];
	}

	inline void sinCosBounds(
		float fIntegral, float *leftSin, float *rightSin, float *leftCos, float *rightCos
	) {
		return sinCosBounds((int32_t)fIntegral, leftSin, rightSin, leftCos, rightCos);
	}

	inline void sinCosBounds(
		float_4 &fIntegral, float_4 *leftSin, float_4 *rightSin, float_4 *leftCos, float_4 *rightCos
	) {
		int32_t indices[4];
		simd::int32_4(fIntegral).store(indices);

		float leftSinBuf[4];
		float rightSinBuf[4];
		float leftCosBuf[4];
		float rightCosBuf[4];

		sinCosBounds(indices[0], leftSinBuf, rightSinBuf, leftCosBuf, rightCosBuf);
		sinCosBounds(indices[1], leftSinBuf + 1, rightSinBuf + 1, leftCosBuf + 1, rightCosBuf + 1);
		sinCosBounds(indices[2], leftSinBuf + 2, rightSinBuf + 2, leftCosBuf + 2, rightCosBuf + 2);
		sinCosBounds(indices[3], leftSinBuf + 3, rightSinBuf + 3, leftCosBuf + 3, rightCosBuf + 3);

		leftSin->v = _mm_loadu_ps(leftSinBuf);
		rightSin->v = _mm_loadu_ps(rightSinBuf);
		leftCos->v = _mm_loadu_ps(leftCosBuf);
		rightCos->v = _mm_loadu_ps(rightCosBuf);
	}

	template<typename T>
	void sinCos2Pi(T x, T *sinOut, T *cosOut) {
		T fIdx = x * SIZE;
		T fFloor = simd::floor(fIdx);
		fIdx -= fFloor;
		
		T leftSin, rightSin, leftCos, rightCos;
		sinCosBounds(fFloor, &leftSin, &rightSin, &leftCos, &rightCos);

		T leftSlopeMult = FRAC * fIdx;
		T rightSlopeMult = FRAC - leftSlopeMult;
		*sinOut = simd::crossfade(leftSin + leftCos * leftSlopeMult, rightSin - rightCos * rightSlopeMult, fIdx);
		*cosOut = simd::crossfade(leftCos - leftSin * leftSlopeMult, rightCos + rightSin * rightSlopeMult, fIdx);
	}
};