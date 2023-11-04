#pragma once

#include <functional>
#include "rack.hpp"

using namespace rack;

// See: https://ccrma.stanford.edu/~jatin/Notebooks/adaa.html
template<typename T, class S>
struct WaveshapingADAA1 {
    T lastInput;
    T lastAntiderivative;

    WaveshapingADAA1(T lastInputInit = T(0)): lastInput(lastInputInit) {
        this->lastAntiderivative = getAD(lastInputInit);
    }

    T process(T input);
    static T getAD(T input) { return S::antiderivative(input); };
};

template<typename T>
struct QuadraticDistortionADAA1 : WaveshapingADAA1<T, QuadraticDistortionADAA1<T>> {
    static T antiderivative(T input) {
        T abs = simd::abs(input);
		T sign = simd::sgn(input);
        T input2 = input * input;
        T input3 = input2 * input;
        T out = simd::ifelse(abs >= T(2), T(1.5) * sign * input, input2 * T(0.5));
        return simd::ifelse(abs <= T(1), out, input2 * T(2) - sign * (T(0.33333333333) * input3 + T(2.5) * input));
    }
};