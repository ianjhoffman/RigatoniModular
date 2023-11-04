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
struct QuadraticDistortionADAA1 : WaveshapingADAA1<T, QuadraticDistortionADAA1<T>> {
    /*
     * This quadratic distortion is implemented as:
     *
     *                       sgn(x)*1.5  , |x| >= 2
     *                                x  , |x| <= 1
     *                -1.5 + (-x - 2)^2  , 1 < x < 2
     *                  1.5 - (x - 2)^2  , -2 < x < -1
     *
     * so the antiderivative is:
     * 
     *                1.5x*sgn(x) - 4/3  , |x| >= 2
     *                          (x^2)/2  , |x| <= 1
     *      (x^3)/3 + 2x^2 + 2.5x + 4/3  , x < 0
     *     -(x^3)/3 + 2x^2 - 2.5x + 4/3  , x >= 0
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
        T out = simd::ifelse(abs >= T(2), T(1.5) * sign * input - T(1.33333333333), input2 * T(0.5));
        return simd::ifelse(
            abs <= T(1),
            out,
            T(1.33333333333) + (T(2) * input2) - sign * (T(0.33333333333) * input3 + T(2.5) * input)
        );
    }
};