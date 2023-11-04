#include <algorithm>

#include "EuclideanPatternGenerator.hpp"

// https://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c
int numRelevantBits(uint64_t n) {
    if (n == 0) return 1;
    #define S(k) if (n >= (UINT64_C(1) << k)) { i += k; n >>= k; }
    int i = -(n == 0); S(32); S(16); S(8); S(4); S(2); S(1); return i + 1;
    #undef S
}

inline uint64_t concatenateBitmasks(uint64_t a, uint64_t b) {
    return (a << numRelevantBits(b)) | b;
}

template<uint8_t MAX_LENGTH>
inline uint64_t EuclideanPatternGenerator<MAX_LENGTH>::getShiftedPattern(uint8_t length, uint8_t density, uint8_t shift) {
    auto pattern = this->patternTable[length - 1][density];
    return((pattern >> shift) & this->lengthMasks[length - 1]) | (pattern << (length - shift));
}

// Implementation based on https://medium.com/code-music-noise/euclidean-rhythms-391d879494df
template<uint8_t MAX_LENGTH>
uint64_t EuclideanPatternGenerator<MAX_LENGTH>::calculateEuclideanBitmask(uint8_t length, uint8_t density) {
    if (length == density) return 0xffffffffffffffff << (64 - length);

    std::vector<uint64_t> ons(density, 1);
    std::vector<uint64_t> offs(length - density, 0);
    while (offs.size() > 1) {
        int numCombinations = std::min(ons.size(), offs.size());
        std::vector<uint64_t> combined{};
        for (int i = 0; i < numCombinations; i++) {
            combined.push_back(concatenateBitmasks(ons[i], offs[i]));
        }

        if (ons.size() > offs.size()) {
            // Remaining ons become offs
            offs = std::vector<uint64_t>(ons.begin() + offs.size(), ons.end());
        } else if (ons.size() < offs.size()) {
            // Remaining offs stay offs
            offs = std::vector<uint64_t>(offs.begin() + ons.size(), offs.end());
        } else {
            offs.clear();
        }

        ons = combined;
    }

    uint64_t accum = 0;
    for (auto &&onPattern : ons) {
        accum = concatenateBitmasks(accum, onPattern);
    }
    for (auto &&offPattern : offs) {
        accum = concatenateBitmasks(accum, offPattern);
    }

    // Make first step of pattern the most significant bit
    return accum << (64 - length);
}