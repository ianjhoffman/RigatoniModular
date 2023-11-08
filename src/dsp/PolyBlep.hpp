#pragma once

#include "rack.hpp"

using namespace rack;

/*
 * See:
 *  - https://www.kvraudio.com/forum/viewtopic.php?t=398553
 *  - https://github.com/pichenettes/stmlib/blob/master/dsp/polyblep.h
 */
template<typename T>
struct PolyBlep {
    T thisSample;
    T nextSample;

    PolyBlep(): thisSample(0), nextSample(0) {}

    inline static float thisBlepSample(float t) {
        const float t2 = t * t;
        return t * t2 - 0.5f * t2 * t2;
    }

    inline static float nextBlepSample(float t) {
        return -thisBlepSample(1.f - t);
    }

    inline static float nextIntegratedBlepSample(float t) {    
        constexpr float COEFF_4 = 1.f / 3.f;
        constexpr float COEFF_5 = 1.f / 8.f;

        const float t2 = t * t;
        const float t4 = t2 * t2;
        return t4 * COEFF_4 - t * t4 * COEFF_5;
    }

    inline static float thisIntegratedBlepSample(float t) {
        return nextIntegratedBlepSample(1.f - t);
    }

    void insertDiscontinuity(float t, T magnitude) {
        this->thisSample += thisBlepSample(t) * magnitude;
        this->nextSample += nextBlepSample(t) * magnitude;
    }

    void insert1stDerivativeDiscontinuity(float t, T magnitude) {
        this->thisSample += thisIntegratedBlepSample(t) * magnitude;
        this->nextSample += nextIntegratedBlepSample(t) * magnitude;
    }

    void registerNextSampleValue(T next) {
        this->nextSample += next;
    }

    void step(T *out) {
        *out = this->thisSample;
        this->thisSample = this->nextSample;
        this->nextSample = T(0);
    }
};