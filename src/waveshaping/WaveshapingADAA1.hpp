#pragma once

#include <functional>
#include "rack.hpp"

using namespace rack;

// See: https://ccrma.stanford.edu/~jatin/Notebooks/adaa.html
template<typename T, class S>
struct WaveshapingADAA1 {
    const float DIFF_LIMIT = 1e-5;

    T lastInput;
    T lastAntiderivative;

    WaveshapingADAA1(T lastInputInit = T(0)): lastInput(lastInputInit) {
        this->lastAntiderivative = S::transformAD1(lastInputInit);
    }

    T process(T input) {
        T diff = input - this->lastInput;
        T fallback = simd::abs(diff) < T(DIFF_LIMIT);
        T currAntiderivative = S::transformAD1(input);
        T ret = simd::ifelse(
            fallback,
            S::transform(T(0.5) * (input + this->lastInput)),
            (currAntiderivative - this->lastAntiderivative) * simd::rcp(diff)
        );
        
        this->lastInput = input;
        this->lastAntiderivative = currAntiderivative;
        return ret;
    }
};

template<typename T>
struct HardClipDistortionADAA1 : WaveshapingADAA1<T, HardClipDistortionADAA1<T>> {
    /*
     * The hard clipping transform is implemented as:
     *
     *                    x  , |x| <= 1
     *               sgn(x)  , otherwise
     *
     * and its first antiderivative is:
     * 
     *              (x^2)/2  , |x| <= 1
     *     sgn(x) * x - 0.5  , otherwise
     */
    static T transform(T input) {
        T abs = simd::abs(input);
		T sign = simd::sgn(input);
        return simd::ifelse(abs <= T(1), input, sign);
    }

    static T transformAD1(T input) {
        T abs = simd::abs(input);
		T sign = simd::sgn(input);
        return simd::ifelse(abs <= T(1), input * input * T(0.5), input * sign - T(0.5));
    }
};

template<typename T>
struct QuadraticDistortionADAA1 : WaveshapingADAA1<T, QuadraticDistortionADAA1<T>> {
    /*
     * This quadratic distortion transform is implemented as:
     *
     *                               sgn(x)*1.5  , |x| >= 2
     *                                        x  , |x| <= 1
     *         sgn(x) * (1.5 + 0.5*(|x| - 2)^2)  , otherwise
     *
     * and its first antiderivative is:
     * 
     *                        1.5x*sgn(x) - 1/6  , |x| >= 2
     *                                  (x^2)/2  , |x| <= 1
     *     1/6 + x^2 - sgn(x) * ((x^3)/6 + x/2)  , otherwise
     */
    static T transform(T input) {
		T abs = simd::abs(input);
		T sign = simd::sgn(input);
		T absShifted = abs - 2.f;
        T out = simd::ifelse(abs >= T(2), T(1.5) * sign, input);
        return simd::ifelse(abs <= T(1), out, sign * (T(1.5) - T(0.5) * absShifted * absShifted));
    }

    static T transformAD1(T input) {
        T abs = simd::abs(input);
		T sign = simd::sgn(input);
        T input2 = input * input;
        T input3 = input2 * input;
        T out = simd::ifelse(abs >= T(2), T(1.5) * sign * input - T(0.1666666666), input2 * T(0.5));
        return simd::ifelse(
            abs <= T(1),
            out,
            T(0.1666666666) + input2 - sign * (T(0.1666666666) * input3 + T(0.5) * input)
        );
    }
};