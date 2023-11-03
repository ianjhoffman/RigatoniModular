#pragma once

template<typename T>
struct DistortionADAA1 {
    T lastInput;
    T lastAntiderivative;

    DistortionADAA1(T lastInputInit = T(0)) {
        this->lastInput = lastInputInit;
        this->lastAntiderivative = this->antiderivative(lastInputInit);
    }

    T process(T input);

    // The antiderivative of the distortion function
    // to compute using ADAA1. It's helpful to add a
    // comment with the original function for context.
    // See: https://ccrma.stanford.edu/~jatin/Notebooks/adaa.html
    virtual T antiderivative(T input);
};

template<typename T>
struct QuadraticDistortionADAA1 : DistortionADAA1<T> {
    T antiderivative(T input) override;
};