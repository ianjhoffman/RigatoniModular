#pragma once

#include "rack.hpp"

using namespace rack;

// See https://www.native-instruments.com/fileadmin/ni_media/downloads/pdf/SineSync.pdf
template<int Z, int O, typename T>
struct HighOrderLinearBlep {
    static constexpr int TABLE_SIZE = 2 * Z * O + 1;
    static constexpr int BUF_SIZE = 2 * Z;
    static constexpr int MID_IDX = (TABLE_SIZE - 1) / 2;
    static constexpr double X_SCALE = (double)1 / (double)O * 2;

    float resid0[TABLE_SIZE];
    float resid1[TABLE_SIZE];
    float resid2[TABLE_SIZE];
    float resid3[TABLE_SIZE];

    T buf[BUF_SIZE] = {};
    int readPos = 0;
    int writePos = Z;

    inline static double blackmanWindow(double piX) {
        return 0.58 + 0.5 * std::cos(piX / Z) - 0.08 * std::cos(2.0 * piX / Z);
    }

    inline void populateResiduals() {
        // Calculate 0th order blep, starting by filling table with sinc
        double calcResid[TABLE_SIZE];
        calcResid[MID_IDX] = 1;
        for (int i = 1; i <= MID_IDX; i++) {
            double piX = (M_PI * (double)i) * X_SCALE;
            double windowedSincVal = blackmanWindow(piX) * std::sin(piX) / piX;

            // Sinc is symmetrical so put this on both sides of 0
            calcResid[MID_IDX + i] = windowedSincVal;
            calcResid[MID_IDX - i] = windowedSincVal;
        }

        // Integrate sinc
        for (int i = 1; i < TABLE_SIZE; i++) {
            calcResid[i] += calcResid[i - 1];
        }

        // Normalize, subtract trivial step to get ∆h_0
        float norm = 1.f / calcResid[TABLE_SIZE - 1];
        for (int i = 0; i < TABLE_SIZE; i++) {
            calcResid[i] = (calcResid[i] * norm) - 0.5;
            resid0[i] = calcResid[i] + ((i >= MID_IDX) ? -0.5 : 0.5);
        }

        // Calculate ∆h_1, ∆h_2, ∆h_3 based on ∆h_0
        constexpr double oneThird = (double)1.f / (double)3.f;
        for (int i = 0; i < TABLE_SIZE; i++) {
            double piX = (M_PI * (double)(i - MID_IDX)) * X_SCALE;

            // ∆h_1
            calcResid[i] = piX * calcResid[i] + M_1_PI * std::cos(piX);
            resid1[i] = calcResid[i];

            // ∆h_2
            calcResid[i] = 0.5 * (piX * calcResid[i] + M_1_PI * std::sin(piX));
            resid2[i] = calcResid[i];

            // ∆h_3
            calcResid[i] = oneThird * (piX * calcResid[i] - M_1_PI * std::cos(piX));
            resid3[i] = calcResid[i];

            DEBUG("%f, %f, %f, %f", resid0[i], resid1[i], resid2[i], resid3[i]);
        }
    }

    HighOrderLinearBlep() {
        populateResiduals();
    }
};