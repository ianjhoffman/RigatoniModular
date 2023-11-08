#pragma once

#include "rack.hpp"

using namespace rack;

template<int Z, int O, typename T>
struct HighOrderLinearBlep {
    static const int TABLE_SIZE = 2 * Z * O + 1;
    static const int BUF_SIZE = 2 * Z;

    float resid0[TABLE_SIZE];
    float resid1[TABLE_SIZE];
    float resid2[TABLE_SIZE];
    float resid3[TABLE_SIZE];

    T buf[BUF_SIZE];
    int readPos;
    int writePos;

    HighOrderLinearBlep() {
        // TODO
    }
};