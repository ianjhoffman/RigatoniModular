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

    inline void populateResiduals() {
        // Calculate 0th order blep, starting by filling table with sinc
        double calcResid[TABLE_SIZE];
        calcResid[MID_IDX] = 1;
        for (int i = 1; i <= MID_IDX; i++) {
            double piX = M_PI * (double)i * X_SCALE;
            double windowedSincVal = blackmanWindow(piX) * std::sin(piX) / piX;

            // Sinc is symmetrical so put this on both sides of 0
            calcResid[MID_IDX + i] = windowedSincVal;
            calcResid[MID_IDX - i] = windowedSincVal;
        }

        // Integrate sinc
        calcResid[MID_IDX] = 0;
        for (int i = 1; i <= MID_IDX; i++) {
            double integral = calcResid[MID_IDX + i - 1] + calcResid[MID_IDX + i] * X_SCALE;
            calcResid[MID_IDX + i] = integral;
            calcResid[MID_IDX - i] = -integral;
        }

        // Normalize, subtract trivial step to get ∆h_0
        double norm = 1.0 / (calcResid[TABLE_SIZE - 1] - calcResid[0]);
        for (int i = 0; i < TABLE_SIZE; i++) {
            calcResid[i] = calcResid[i] * norm + ((i >= MID_IDX) ? -0.5 : 0.5);
            resid0[i] = calcResid[i];
        }

        // Calculate ∆h_1, ∆h_2, ∆h_3 based on ∆h_0
        constexpr double oneThird = (double)1 / (double)3;
        for (int i = 0; i < TABLE_SIZE; i++) {
            double piX = (M_PI * (double)(i - MID_IDX)) * X_SCALE;

            // ∆h_1
            calcResid[i] = piX * calcResid[i] + M_1_PI * std::cos(piX);
            calcResid[i] *= blackmanWindow(piX);
            resid1[i] = calcResid[i];

            // ∆h_2
            calcResid[i] = 0.5 * (piX * calcResid[i] + M_1_PI * std::sin(piX));
            calcResid[i] *= blackmanWindow(piX);
            resid2[i] = calcResid[i];

            // ∆h_3
            calcResid[i] = oneThird * (piX * calcResid[i] - M_1_PI * std::cos(piX));
            calcResid[i] *= blackmanWindow(piX);
            resid3[i] = calcResid[i];
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

			buf[bufIndex] -= order0 * interpolateResidual(resid0, blepIndex, true);
            buf[bufIndex] -= order1 * interpolateResidual(resid1, blepIndex);
            buf[bufIndex] -= order2 * interpolateResidual(resid2, blepIndex);
            buf[bufIndex] -= order3 * interpolateResidual(resid3, blepIndex);

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