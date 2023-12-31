#pragma once

// See https://www.native-instruments.com/fileadmin/ni_media/downloads/pdf/SineSync.pdf
template<int Z, int O, typename T>
struct HighOrderLinearBlep {
    static constexpr int TABLE_SIZE = 2 * Z * O + 1;
    static constexpr int BUF_SIZE = 2 * Z;
    static constexpr int MID_IDX = (TABLE_SIZE - 1) / 2;
    static constexpr double X_SCALE = (double)1 / (double)(O * 2);

    float resid0[TABLE_SIZE];
    float resid1[TABLE_SIZE];
    float resid2[TABLE_SIZE];
    float resid3[TABLE_SIZE];

    T buf[BUF_SIZE] = {};
    int readPos = 0;
    int writePos = Z;

    inline static double blackmanWindow(double piX) {
        return 0.58 + 0.5 * std::cos(2.0 * piX / Z) - 0.08 * std::cos(4.0 * piX / Z);
    }

    inline static double trapezoid(double width, double h1, double h2) {
        return width * 0.5 * (h1 + h2);
    }

    inline void populateResiduals() {
        double piXTab[TABLE_SIZE];
        for (int i = 0; i < TABLE_SIZE; i++) {
            piXTab[i] = (M_PI * (double)(i - MID_IDX)) * X_SCALE;
        }

        // Calculate 0th order blep, starting by filling table with sinc
        double calcSinc[TABLE_SIZE];
        double calcResid[TABLE_SIZE];
        calcSinc[MID_IDX] = 1;
        for (int i = 1; i <= MID_IDX; i++) {
            double piX = piXTab[MID_IDX + i];
            double sincVal = std::sin(piX) / piX;

            // Sinc is symmetrical so put this on both sides of 0
            calcSinc[MID_IDX + i] = sincVal;
            calcSinc[MID_IDX - i] = sincVal;
        }

        // Integrate sinc to get Si(x)
        calcResid[MID_IDX] = 0;
        for (int i = 1; i <= MID_IDX; i++) {
            int j = MID_IDX + i;
            double slice = trapezoid(X_SCALE, calcSinc[j - 1], calcSinc[j]);
            double integral = calcResid[j - 1] + slice;
            calcResid[j] = integral;
            calcResid[MID_IDX - i] = -integral; // invert for x < 0
        }

        // Calculate all 4 residuals
        for (int i = 0; i < TABLE_SIZE; i++) {
            double piX = piXTab[i];
            double blackmanMult = blackmanWindow(piX);

            // Subtract trivial step from sine integral to get ∆h_0
            calcResid[i] = calcResid[i] + ((i >= MID_IDX) ? -0.5 : 0.5);
            resid0[i] = blackmanMult * calcResid[i];

            // Calculate ∆h_1 from ∆h_0
            calcResid[i] = piX * calcResid[i] + M_1_PI * std::cos(piX);
            resid1[i] = blackmanMult * calcResid[i];

            // Calculate ∆h_2 from ∆h_1
            calcResid[i] = (piX * calcResid[i] + M_1_PI * std::sin(piX)) * 0.5;
            resid2[i] = blackmanMult * calcResid[i];

            // Calculate ∆h_3 from ∆h_2
            calcResid[i] = (piX * calcResid[i] - M_1_PI * std::cos(piX)) * 0.3333333333333333;
            resid3[i] = blackmanMult * calcResid[i];
        }
    }

    HighOrderLinearBlep() {
        populateResiduals();
    }

    /*
     * Similar to simple linear interpolation of an array
     * but with special rules to not interpolate across any
     * discontinuities at the midpoint of our residual tables.
     */
    float interpolateResidual(float *resid, float fracIndex, bool flipAtCenter = false) {
        int idx = fracIndex;
        float fade = fracIndex - idx;

        float a = resid[idx];
        float b = ((idx == MID_IDX - 1) && flipAtCenter) ? -resid[idx + 1] : resid[idx + 1];
        return a + (b - a) * fade;
    }

    /*
     * Insert discontinuities in the 0th, 1st, 2nd, and 3rd derivatives of the signal to
     * be smoothed by BLEP. The discontinuity is centered -1 < t <= 0 samples in the past.
     */
    void insertDiscontinuities(float t, const T &order0, const T &order1, const T &order2, const T &order3) {
        if (!(0 <= t && t < 1)) return;

        int bufIndex = readPos;
        for (int i = 0; i < BUF_SIZE; i++) {
            float blepIndex = ((float)i + t) * O;

            buf[bufIndex] += order0 * interpolateResidual(resid0, blepIndex, true);
            buf[bufIndex] += order1 * interpolateResidual(resid1, blepIndex);
            buf[bufIndex] += order2 * interpolateResidual(resid2, blepIndex);
            buf[bufIndex] += order3 * interpolateResidual(resid3, blepIndex);

            bufIndex = (bufIndex + 1 == BUF_SIZE) ? 0 : bufIndex + 1;
        }
    }

    void insertDiscontinuities(float t, const T &order0, const T &order1, const T &order2) {
        const T NO_DISC = T(0);
        insertDiscontinuities(t, order0, order1, order2, NO_DISC);
    }

    void insertDiscontinuities(float t, const T &order0, const T &order1) {
        const T NO_DISC = T(0);
        insertDiscontinuities(t, order0, order1, NO_DISC, NO_DISC);
    }

    void insertDiscontinuities(float t, const T &order0) {
        const T NO_DISC = T(0);
        insertDiscontinuities(t, order0, NO_DISC, NO_DISC, NO_DISC);
    }

    /*
     * Write the current uncorrected sample to a time slightly in the future and
     * return the current corrected sample to the caller in the provided reference.
     */
    void processSample(T &currSample, T &out) {
        out = buf[readPos];
        buf[readPos] = T(0);
        buf[writePos] += currSample;
        readPos = (readPos + 1) % BUF_SIZE;
        writePos = (writePos + 1) % BUF_SIZE;
    }
};