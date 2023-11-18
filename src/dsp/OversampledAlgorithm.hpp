#include "rack.hpp"

#include <algorithm>
#include <array>

using namespace rack;

enum ParameterInterpolationType {
    LINEAR
    // TODO - cubic
};

template<int NUM_PARAMS, typename T = float>
struct ParameterInterpolator {
    T paramBuffer[NUM_PARAMS * 2];
    int prevOffset = 0;
    int currOffset = NUM_PARAMS;
    ParameterInterpolationType interpolationType;

    ParameterInterpolator(ParameterInterpolationType type = ParameterInterpolationType::LINEAR) {
        this->interpolationType = type;
        std::fill(paramBuffer, paramBuffer + NUM_PARAMS * 2, T(0));
    }

    void update(std::array<T, NUM_PARAMS> &newParams) {
        std::swap(this->prevOffset, this->currOffset);
        std::copy(newParams.begin(), newParams.end(), paramBuffer + currOffset);
    }

    void interpolate(std::array<T, NUM_PARAMS> &dest, float fade) {
        for (int i = 0; i < NUM_PARAMS; i++) {
            switch (this->interpolationType) {
                case ParameterInterpolationType::LINEAR:
                    auto prev = paramBuffer[prevOffset + i];
                    auto curr = paramBuffer[currOffset + i];
                    dest[i] = prev + (curr - prev) * T(fade);
            }
        }
    }
};

template<int OVERSAMPLE, int QUALITY, int NUM_OUTS, int NUM_PARAMS, typename TParams = float, typename TOut = float>
struct OversampledAlgorithm {
    ParameterInterpolator<NUM_PARAMS, TParams> paramInterp;
    dsp::Decimator<OVERSAMPLE, QUALITY, TOut> decimators[NUM_OUTS];
    float step;
    bool oversamplingEnabled{true};

    OversampledAlgorithm(ParameterInterpolator<NUM_PARAMS, TParams> interp): paramInterp(interp) {
        step = 1.f / (float)OVERSAMPLE;
    }

    int getMultiplier() {
        return oversamplingEnabled ? OVERSAMPLE : 1;
    }

    float getDivisor() {
        return oversamplingEnabled ? this->step : 1.f;
    }

    std::array<TOut, NUM_OUTS> process(const Module::ProcessArgs& args, std::array<TParams, NUM_PARAMS> &params) {
        if (!oversamplingEnabled) return this->processFrame(args, params);

        // Modify args to account for oversampled sample rate
        Module::ProcessArgs processArgs = {
            .sampleRate = args.sampleRate * OVERSAMPLE,
            .sampleTime = args.sampleTime * this->step,
            .frame = args.frame
        };

        this->paramInterp.update(params);
        float progress = this->step;
        TOut decimatorInputs[NUM_OUTS][OVERSAMPLE];
        std::array<TParams, NUM_PARAMS> stepParams;
        for (int i = 0; i < OVERSAMPLE; i++) {
            this->paramInterp.interpolate(stepParams, progress);
            auto frameOuts = this->processFrame(processArgs, stepParams);
            for (int j = 0; j < NUM_OUTS; j++) {
                decimatorInputs[j][i] = frameOuts[j];
            }

            progress += step;
        }

        std::array<TOut, NUM_OUTS> out;
        for (int i = 0; i < NUM_OUTS; i++) {
            out[i] = this->decimators[i].process(decimatorInputs[i]);
        }

        return out;
    }

    virtual std::array<TOut, NUM_OUTS> processFrame(const Module::ProcessArgs& args, std::array<TParams, NUM_PARAMS> &params) = 0;
};