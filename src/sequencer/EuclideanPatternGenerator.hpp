#pragma once

#include <utility>

struct EuclideanPatternGenerator {
    // All Euclidean pattern bitmasks for each length up to 64; inner index is density
	std::array<std::array<uint64_t, 64>, 64> patternTable;

    // Masks for shifting sequences of each length
	std::array<uint64_t, 64> lengthMasks;

    EuclideanPatternGenerator() {
        for (int length = 1; length <= 64; length++) {
			for (int density = 1; density <= length; density++) {
				this->patternTable[length - 1][density - 1] = calculateEuclideanBitmask(length, density);
			}
			this->lengthMasks[length - 1] = 0xffffffffffffffff << (64 - length);
		}
    }

    static uint64_t calculateEuclideanBitmask(uint8_t length, uint8_t density);
    uint64_t getShiftedPattern(uint8_t length, uint8_t density, uint8_t shift);
};