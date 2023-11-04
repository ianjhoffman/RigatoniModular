#include "rack.hpp"
#include "WaveshapingADAA1.hpp"

using namespace rack;

constexpr float DIFF_LIMIT = 1e-5;

template<typename T, class S>
T WaveshapingADAA1<T, S>::process(T input) {
    T diff = input - this->lastInput;
    T fallback = simd::abs(diff) < DIFF_LIMIT;
    T currAntiderivative = getAD(input);
    T ret = simd::ifelse(
        fallback,
        T(0.5) * (input + this->lastInput),
        (currAntiderivative - this->lastAntiderivative) * simd::rcp(diff)
    );
    
    this->lastInput = input;
    this->lastAntiderivative = currAntiderivative;
    return ret;
}

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
template<typename T>
T quadraticDistortionAD(T input) {
    	T abs = simd::abs(input);
		T sign = simd::sgn(input);
        T input2 = input * input;
        T input3 = input2 * input;
        T out = simd::ifelse(abs >= T(2), T(1.5) * sign * input, input2 * T(0.5));
        return simd::ifelse(abs <= T(1), out, input2 * T(2) - sign * (T(0.33333333333) * input3 + T(2.5) * input));
}

// TODO: consider implementing ADAA1 versions of these distortions:

/*

// Soft clipping (quadratic) drive with a small amount of wave folding at the extreme
static float_4 drive2(float_4 in, float drive) {
    // Don't hit the folder as hard, up to 12 o'clock on the drive knob should almost never clip
    in *= drive * .5f;
    auto overThresh = simd::abs(in) - .9f;
    auto underThreshMask = overThresh < 0.f;
    auto subAmt = overThresh * overThresh * drive * simd::sgn(in);
    return simd::ifelse(underThreshMask, in, in - subAmt);
}

// Cubic soft clipping with no foldback
static float_4 drive3(float_4 in, float drive) {
    constexpr float BASE_SLOPE = 0.078536469359214f; // 3/32 * (4 - sqrt(10))
    constexpr float CURVE_ADD = -0.988211768802619f; // 1/4 - (5/8 * sqrt(5/2))
    constexpr float X_SCALE = 8.f;
    constexpr float X_LIMIT = 16.f; // Slope of soft clipping function is 0 at x = 16
    constexpr float Y_AT_LIMIT = 1.011788231197381f; // func(16) = 2 - (5/8 * sqrt(5/2))

    auto x = in * X_SCALE * drive;
    auto abs = simd::abs(x);
    auto sign = simd::sgn(x);
    auto underThreshMask = abs <= 10.f;
    auto addAmt = sign * (-0.0625f * simd::pow(abs, float_4(1.5f)) + CURVE_ADD);
    auto out = simd::ifelse(underThreshMask, BASE_SLOPE * x, 0.375f * x + addAmt);
    return simd::ifelse(abs >= X_LIMIT, sign * Y_AT_LIMIT, out);
}

*/