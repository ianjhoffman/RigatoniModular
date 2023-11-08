#include "plugin.hpp"
#include "sequencer/EuclideanPatternGenerator.hpp"
#include "math/Utility.hpp"
#include "dsp/OversampledAlgorithm.hpp"
#include "dsp/PolyBlep.hpp"
#include "dsp/WaveshapingADAA1.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <utility>

using simd::float_4;

enum ContinuousStrideMode {
	OFF,
	SYNC,
	FREE
};

const uint64_t AMP_MASK = 0x000000000000000f;
const int AMP_SHIFT = 60;
const float LFO_MULTIPLIER = .05f;
const float VCO_MULTIPLIER = 20.f;
const float LIN_FM_FACTOR = 5.f;
const float MAX_FREQ = 10000.f;
constexpr float VCO_FINE_TUNE_OCTAVES = 7.f / 12.f; // a fifth

// Swap order of every 4-bit chunk
inline static uint64_t flipNibbleEndian(uint64_t a) {
	uint64_t out = (a & 0x8888888888888888) >> 3; // 3 to 0
	out |= (a & 0x1111111111111111) << 3; // 0 to 3
	out |= (a & 0x4444444444444444) >> 1; // 2 to 1
	out |= (a & 0x2222222222222222) << 1; // 1 to 2
	return out;
}

float scaleLength(float knobValue) {
	return clamp(dsp::exp2_taylor5(knobValue) * 16.8f - 3.2f, 1.f, 64.f);
}

float scaleStride(float knobValue) {
	return clamp(dsp::exp2_taylor5(knobValue) * 1.333333333333f - 1.333333333333f, 0.f, 4.f);
}

void shapeAmplitudes(
	std::array<float_4, 16> &amplitudes,
	uint64_t harmonicMask,
	int numBlocks,
	float length,
	int tilt, // 0-2 (lowpass, bandpass, highpass)
	float pivotHarm,
	float intensity
) {
	constexpr float SLOPE_SCALE = 0.01f;
	constexpr float SLOPES[3][2] = {
		{                0.f, 10.f * SLOPE_SCALE},
		{-15.f * SLOPE_SCALE, 15.f * SLOPE_SCALE},
		{-10.f * SLOPE_SCALE, 0.f},
	};

	float_4 belowSlope = SLOPES[tilt][0] * intensity;
	float_4 aboveSlope = SLOPES[tilt][1] * intensity;
	float_4 pivotBase = crossfade(0.f, (tilt == 1) ? .25f : 0.f, intensity);
	float_4 pivotDiffs = pivotHarm - float_4(0.f, 1.f, 2.f, 3.f);
	auto shiftAmt = AMP_SHIFT;
	harmonicMask = flipNibbleEndian(harmonicMask);
	for (int i = 0; i < numBlocks; i++) {
		float_4 slopes = simd::ifelse(pivotDiffs > 0, belowSlope, aboveSlope);
		auto addMask = simd::movemaskInverse<float_4>((harmonicMask >> shiftAmt) & AMP_MASK);
		auto toAdd = simd::ifelse(addMask, pivotBase + slopes * pivotDiffs, 0.f);
		amplitudes[i] = clamp(amplitudes[i] + toAdd, 0.f, 1.5f); // No negative amplitudes

		pivotDiffs -= 4.f;
		shiftAmt -= 4;
	}
}

inline float calculateFrequencyHz(float expCv, float linCv, float multiplier) {
	// Clamp to max frequency, allow negative frequency for thru-zero linear FM
	return std::min(multiplier * (dsp::exp2_taylor5(expCv) + linCv), MAX_FREQ);
}

struct LoomAlgorithm : OversampledAlgorithm<2, 10, 1, 3, float_4, float_4> {
	// Anti-aliased drive processor
	QuadraticDistortionADAA1<float_4> driveProcessor{};

	// Euclidean pattern generator/storage
	EuclideanPatternGenerator patternGenerator{};

	// Masks for which harmonics go to which output, based on output mode
	std::array<std::array<float_4, 16>, 3> harmonicSplitMasks;

	// Precalculated amplitude multipliers per harmonic to avoid lots of divisions during runtime
	std::array<float_4, 16> baseAmplitudes;

	// Synthesis parameters
	std::array<float_4, 16> phaseAccumulators;
	PolyBlep<float_4> syncBlep;
	dsp::MinBlepGenerator<16, 16, float> squareBlep;
	float lastSyncValue{0.f};
	float freqMultiplier{VCO_MULTIPLIER};

	// Convenient place to store phases during calculations
	std::array<float_4, 16> sinPhaseWithoutSync;
	std::array<float_4, 16> cosPhaseWithoutSync;
	std::array<float_4, 16> sinPhaseWithSync;
	std::array<float_4, 16> cosPhaseWithSync;

	// Discrete switched parameters that aren't interpolated
	bool splitMode{false};
	bool boostFund{true};
	int tilt{0};
	ContinuousStrideMode continuousStrideMode{ContinuousStrideMode::OFF};

	// Readable parameters for lights
	bool oscLight;
	bool strideLight;
	float_4 ampLights[8];

	// Additional configurable parameters for anti-aliasing
	bool doADAA{true};
	bool doBlep{true};

	LoomAlgorithm() : OversampledAlgorithm<2, 10, 1, 3, float_4, float_4>(ParameterInterpolator<3, float_4>()) {
		// Set up everything pre-cached/pre-allocated
		this->populateHarmonicSplitMasks();
		this->calculateBaseAmplitudes();
		this->phaseAccumulators.fill(0.f);		
	}

	// Gets the Euclidean pattern with the supplied parameters, masked and with
	// each nibble flipped for simd::ifelse (which expects the opposite order)
	inline uint64_t getPattern(uint64_t mask, uint8_t length, uint8_t density, uint8_t shift) {
		return flipNibbleEndian(mask & this->patternGenerator.getShiftedPattern(length, density, shift));
	}

	void populateHarmonicSplitMasks() {
		// 0 = Split mode not active
		this->harmonicSplitMasks[0].fill(simd::movemaskInverse<float_4>(0b1111));

		// 1 = Split mode active, odd harmonics
		this->harmonicSplitMasks[1].fill(simd::movemaskInverse<float_4>(0b0101));

		// 2 = Split mode active, even harmonics (still include fundamental)
		this->harmonicSplitMasks[2].fill(simd::movemaskInverse<float_4>(0b1010));
		this->harmonicSplitMasks[2][0] = simd::movemaskInverse<float_4>(0b1011);
	}

	void calculateBaseAmplitudes() {
		for (int i = 0; i < 64; i++) {
			int idx = i >> 2;
			int offset = i & 0b11;
			this->baseAmplitudes[idx][offset] = 1.f / (i + 1);
		}
	}

	inline void setAmplitudesHalf(
		std::array<float_4, 16> &amplitudes,
		int length,
		float lengthFade,
		float density,
		float shift,
		uint64_t harmonicMask,
		int numBlocks
	) {
		int lengthIdx = length - 1;

		// Density patterns and fade amount
		float fDensity = density * lengthIdx;
		int iDensityLow = (int)std::floor(fDensity);
		int iDensityHigh = std::min(iDensityLow + 1, lengthIdx);
		float densityFade = fDensity - iDensityLow;

		// Shift values and fade amount
		float fShift = shift * lengthIdx;
		int iShiftLow = (int)std::floor(fShift);
		int iShiftHigh = iShiftLow + 1;
		float shiftFade = fShift - iShiftLow;

		float_4 fader = {
			(1.f - densityFade) * (1.f - shiftFade), // low density, low shift
			(1.f - densityFade) * shiftFade, // low density, high shift
			densityFade * (1.f - shiftFade), // high density, low shift
			densityFade * shiftFade // high density, high shift
		};
		fader *= lengthFade;

		// Do blending
		auto lowDensLowShift = this->getPattern(harmonicMask, length, iDensityLow, iShiftLow);
		auto lowDensHighShift = this->getPattern(harmonicMask, length, iDensityLow, iShiftHigh);
		auto highDensLowShift = this->getPattern(harmonicMask, length, iDensityHigh, iShiftLow);
		auto highDensHighShift = this->getPattern(harmonicMask, length, iDensityHigh, iShiftHigh);
		auto shiftAmt = AMP_SHIFT;
		for (int i = 0; i < numBlocks; i++) {
			auto lowLow = simd::movemaskInverse<float_4>((lowDensLowShift >> shiftAmt) & AMP_MASK);
			auto lowHigh = simd::movemaskInverse<float_4>((lowDensHighShift >> shiftAmt) & AMP_MASK);
			auto highLow = simd::movemaskInverse<float_4>((highDensLowShift >> shiftAmt) & AMP_MASK);
			auto highHigh = simd::movemaskInverse<float_4>((highDensHighShift >> shiftAmt) & AMP_MASK);

			// Wish I could do a matrix multiplication here
			amplitudes[i] += simd::ifelse(lowLow, fader[0], 0.f);
			amplitudes[i] += simd::ifelse(lowHigh, fader[1], 0.f);
			amplitudes[i] += simd::ifelse(highLow, fader[2], 0.f);
			amplitudes[i] += simd::ifelse(highHigh, fader[3], 0.f);

			shiftAmt -= 4;
		}
	}

	void setAmplitudes(
		std::array<float_4, 16> &amplitudes,
		int harmonicLimit,
		float length,  // 1-64
		float density, // 0-1
		float stride,  // 0-4
		float shift,   // 0-1
		int *numBlocks,
		uint64_t *harmonicMask
	) {
		*harmonicMask = 0xffffffffffffffff << (64 - harmonicLimit);
		int iLengthLow = (int)std::floor(length);
		int iLengthHigh = std::min(iLengthLow + 1, 64);
		*numBlocks = (std::min(harmonicLimit, iLengthHigh) + 0b11) >> 2;
		float lengthFade = length - iLengthLow;

		this->setAmplitudesHalf(amplitudes, iLengthHigh, lengthFade, density, shift, *harmonicMask, *numBlocks);
		this->setAmplitudesHalf(amplitudes, iLengthLow, (1.f - lengthFade), density, shift, *harmonicMask, *numBlocks);
	}

	std::array<float_4, 1> processFrame(const Module::ProcessArgs& args, std::array<float_4, 3> &params) override {
		// Unpack params
		float length = params[0][0];
		float density = params[0][1];
		float shift = params[0][2];
		float stride = params[0][3];

		float expCv = params[1][0];
		float linCv = params[1][1];
		float sinPhaseOffset = params[1][2];
		float syncValue = params[1][3];

		float pivotHarm = params[2][0];
		float intensity = params[2][1];
		float drive = params[2][2];

		// Perform any necessary modifications on unpacked parameters
		length = scaleLength(clamp(length, -2.f, 2.f));
		density = clamp(density);
		shift = clamp(shift);
		stride = scaleStride(clamp(stride, 0.f, 2.f));
		if (this->continuousStrideMode == ContinuousStrideMode::OFF) stride = std::round(stride);
		this->strideLight = std::abs(stride - 1.f) < .01f; // Reasonable epsilon

		sinPhaseOffset = clamp(sinPhaseOffset, -1.f, 1.f);
		sinPhaseOffset -= std::floor(sinPhaseOffset);

		pivotHarm = scaleLength(clamp(pivotHarm, -2.f, 2.f)) - 1.f;
		intensity = clamp(intensity);
		drive = clamp(drive, 0.f, 2.f);

		// Calculate oscillator frequency based on a multitude of factors
		float freq = calculateFrequencyHz(expCv, linCv, this->freqMultiplier);

		// Calculate harmonic limit for anti-aliasing, as well as associated mask for shifted patterns.
		// This functions as simple built-in band-limiting (not perfect because of drive, fm, pm,
		// all of which can introduce extra harmonics). N = ((sample_rate / 2*freq) - 1) / stride
		// Work off of the sample rate pre-oversampling by multiplying by step
		auto freq2Recip = simd::rcp(2.f * std::abs(freq)); // TZFM can make frequency negative
		float harmonicMultipleLimit = args.sampleRate * this->getDivisor() * freq2Recip[0] - 1.f;
		float clampedStride = std::fmax(stride, 0.1f);
		int harmonicLimit = clamp((int)(harmonicMultipleLimit / clampedStride), 1, 64);

		// Actually do partial amplitude/frequency calculations
		std::array<float_4, 16> harmonicAmplitudes{};
		int numBlocks;
		uint64_t harmonicMask;
		this->setAmplitudes(harmonicAmplitudes, harmonicLimit, length, density, stride, shift, &numBlocks, &harmonicMask);

		// Set summed amplitudes for light show
		for (int i = 0; i < 8; i++) {
			int off = i << 1;
			this->ampLights[i] = harmonicAmplitudes[off] + harmonicAmplitudes[off + 1];
		}

		// Shaping section
		shapeAmplitudes(harmonicAmplitudes, harmonicMask, numBlocks, length, tilt, pivotHarm, intensity);

		// Fundamental boosting
		harmonicAmplitudes[0][0] = std::max(harmonicAmplitudes[0][0], this->boostFund ? .5f : 0.f);

		// Calculate sync
		float deltaSync = syncValue - this->lastSyncValue;
		float syncCrossing = -this->lastSyncValue / deltaSync;
		this->lastSyncValue = syncValue;
		bool normalSync = (0.f < syncCrossing) && (syncCrossing <= 1.f) && (syncValue >= 0.f);

		// Calculate phase modulation offsets
		float cosPhaseOffset = sinPhaseOffset + 0.25f;
		cosPhaseOffset -= std::floor(cosPhaseOffset);

		// Additive synthesis params setup
		float phaseInc = freq * args.sampleTime;
		phaseInc -= std::floor(phaseInc);
		float normalSyncPhase = (1.f - syncCrossing) * phaseInc;
		float_4 harmonicFadeLowerBound = harmonicMultipleLimit * .666666666f;
		float_4 harmonicFalloffSlope = simd::rcp(-.333333333333f * harmonicMultipleLimit);

		// Calculate fundamental sync for non-free continuous stride
		float fundAccumBeforeInc = this->phaseAccumulators[0][0];
		float fundPhaseWrapped = fundAccumBeforeInc + phaseInc - 1.f;
		bool fundSync = fundPhaseWrapped >= 0.f;

		// 2 types of sync that can introduce discontinuities
		bool doSync = false;
		float syncPhase = 0.f;
		float minBlepP = 0.f;
		if (normalSync) {
			doSync = true;
			minBlepP = syncCrossing;
			syncPhase = normalSyncPhase;
		} else if (continuousStrideMode != ContinuousStrideMode::FREE && fundSync) {
			doSync = true;
			minBlepP = 1.f - (fundPhaseWrapped / phaseInc);
			syncPhase = fundPhaseWrapped;
		}

		// Sadly for band-limiting with MinBLEP we have to calculate all the partials
		// for each output for both synced and unsynced sine and cosine
		std::array<float_4, 16> &oddHarmSplitMask = this->harmonicSplitMasks[this->splitMode ? 1 : 0];
		std::array<float_4, 16> &evenHarmSplitMask = this->harmonicSplitMasks[this->splitMode ? 2 : 0];
		float_4 out1WithoutSyncSum{}, out1WithSyncSum{}, out2WithoutSyncSum{}, out2WithSyncSum{};
		float_4 out1DerivativeDiscSum{}, out2DerivativeDiscSum{};

		// Trying something new with progressive phase multiple calculations.
		// Calculate the base values for various synthesis parameters for the first 4 harmonics
		// as well as the increments of those values for every block of 4 harmonics. This lets
		// us get rid of a bunch of multiplies inside the synthesis loop.
		float_4 multiples = {1.f, 1.f + stride, 1.f + 2 * stride, 1.f + 3 * stride};
		float_4 basePhaseInc = multiples * phaseInc;
		float_4 out1SyncPhase = multiples * syncPhase;
		float_4 sinPm = multiples * sinPhaseOffset;
		float_4 cosPm = multiples * cosPhaseOffset;

		float_4 phaseIncAdd = phaseInc * 4.f * stride;
		float_4 syncPhaseAdd = syncPhase * 4.f * stride;
		float_4 sinPhaseOffsetAdd = sinPhaseOffset * 4.f * stride;
		float_4 cosPhaseOffsetAdd = sinPhaseOffsetAdd + stride;
		float_4 multipleAdd = 4.f * stride;

		float_4 overallAmplitude, sinWithoutSync, cosWithoutSync;
		for (int i = 0; i < 16; i++) {
			// Always update phase accumulators
			sinPhaseWithoutSync[i] = this->phaseAccumulators[i] + basePhaseInc;
			sinPhaseWithoutSync[i] -= simd::floor(sinPhaseWithoutSync[i]);
			float_4 phaseAccumVal = sinPhaseWithoutSync[i]; // store before doing PM

			if (i < numBlocks) {
				overallAmplitude = harmonicAmplitudes[i] * simd::fmin(simd::rcp(multiples), this->baseAmplitudes[i]);
				// Harmonic falloff between ~.75x and ~1x Nyquist
				overallAmplitude *= simd::clamp(
					1.f + (multiples - harmonicFadeLowerBound) * harmonicFalloffSlope
				);
				
				// Phase modulation
				sinPhaseWithoutSync[i] += sinPm;
				sinPhaseWithoutSync[i] -= simd::floor(sinPhaseWithoutSync[i]);

				cosPhaseWithoutSync[i] = phaseAccumVal + cosPm;
				cosPhaseWithoutSync[i] -= simd::floor(cosPhaseWithoutSync[i]);

				sinWithoutSync = sin2pi_chebyshev(sinPhaseWithoutSync[i]) * overallAmplitude;
				cosWithoutSync = sin2pi_chebyshev(cosPhaseWithoutSync[i]) * overallAmplitude;

				out1WithoutSyncSum += simd::ifelse(oddHarmSplitMask[i], sinWithoutSync, 0.f);
				out2WithoutSyncSum += simd::ifelse(evenHarmSplitMask[i], this->splitMode ? sinWithoutSync : cosWithoutSync, 0.f);
			}

			// It turns out that branching in the loop is significantly faster than branchless
			if (doSync) {
				sinPhaseWithSync[i] = out1SyncPhase;
				sinPhaseWithSync[i] -= simd::floor(sinPhaseWithSync[i]);
				phaseAccumVal = sinPhaseWithSync[i]; // store before doing PM

				if (i < numBlocks) {
					// Phase modulation (sync)
					sinPhaseWithSync[i] += sinPm;
					sinPhaseWithSync[i] -= simd::floor(sinPhaseWithSync[i]);

					cosPhaseWithSync[i] = phaseAccumVal + cosPm;
					cosPhaseWithSync[i] -= simd::floor(cosPhaseWithSync[i]);

					float_4 sinWithSync = sin2pi_chebyshev(sinPhaseWithSync[i]) * overallAmplitude;
					float_4 cosWithSync = sin2pi_chebyshev(cosPhaseWithSync[i]) * overallAmplitude;

					out1WithSyncSum += simd::ifelse(oddHarmSplitMask[i], sinWithSync, 0.f);
					out2WithSyncSum += simd::ifelse(evenHarmSplitMask[i], this->splitMode ? sinWithSync : cosWithSync, 0.f);

					// 1st derivatives
					float_4 sinDiff = sinWithSync - sinWithoutSync;
					float_4 cosDiff = cosWithSync - cosWithoutSync;

					out1DerivativeDiscSum += simd::ifelse(oddHarmSplitMask[i], cosDiff, 0.f);
					out2DerivativeDiscSum += simd::ifelse(evenHarmSplitMask[i], this->splitMode ? cosDiff : -sinDiff, 0.f);
				}
			}

			this->phaseAccumulators[i] = phaseAccumVal;

			// Update progressive phase multiple calculations
			basePhaseInc += phaseIncAdd;
			out1SyncPhase += syncPhaseAdd;
			sinPm += sinPhaseOffsetAdd;
			cosPm += cosPhaseOffsetAdd;
			multiples += multipleAdd;
		}

		// We weren't updating these in the loop so set them here
		if (!doSync) {
			out1WithSyncSum = out1WithoutSyncSum;
			out2WithSyncSum = out2WithoutSyncSum;
		}

		// Collapse the remaining 4 harmonic accumulators for each output
		float_4 mainOutsPacked = {out1WithoutSyncSum[0], out1WithSyncSum[0], out2WithoutSyncSum[0], out2WithSyncSum[0]};
		mainOutsPacked += {out1WithoutSyncSum[1], out1WithSyncSum[1], out2WithoutSyncSum[1], out2WithSyncSum[1]};
		mainOutsPacked += {out1WithoutSyncSum[2], out1WithSyncSum[2], out2WithoutSyncSum[2], out2WithSyncSum[2]};
		mainOutsPacked += {out1WithoutSyncSum[3], out1WithSyncSum[3], out2WithoutSyncSum[3], out2WithSyncSum[3]};
		mainOutsPacked *= (float)this->splitMode + 1.f;

		// Fundamental and square outs are a lot easier
		float_4 fundOutParams = sin2pi_chebyshev<float_4>({
			sinPhaseWithoutSync[0][0], sinPhaseWithSync[0][0],
			cosPhaseWithoutSync[0][0], cosPhaseWithSync[0][0],
		});

		// Oscillator light indicator based on fundamental phase accumulator
		this->oscLight = this->phaseAccumulators[0][0] > 0.5f;

		// Do drive before BLEP, driving BLEP might introduce more aliasing
		mainOutsPacked *= 0.5f * drive;
		mainOutsPacked = this->doADAA ? this->driveProcessor.process(mainOutsPacked) : this->driveProcessor.transform(mainOutsPacked);

		// BLEP
		if (doSync && this->doBlep) {
			float_4 discontinuities = {
				mainOutsPacked[1] - mainOutsPacked[0],
				mainOutsPacked[3] - mainOutsPacked[2],
				fundOutParams[1] - fundOutParams[0],
				0.f,
			};
			float_4 derivativeDiscontinuities = {
				sum_float4(out1DerivativeDiscSum),
				sum_float4(out2DerivativeDiscSum),
				fundOutParams[3] - fundOutParams[2],
				0.f,
			};

			this->syncBlep.insertDiscontinuity(minBlepP, discontinuities);
			this->syncBlep.insert1stDerivativeDiscontinuity(minBlepP, derivativeDiscontinuities);
		}

		this->syncBlep.registerNextSampleValue({mainOutsPacked[1], mainOutsPacked[3], fundOutParams[1], 0.f});
		float_4 outsPacked;
		this->syncBlep.step(&outsPacked);
		return std::array<float_4, 1>{5.f * outsPacked};
	}
};

struct Loom : Module {
	enum ParamId {
		CONTINUOUS_STRIDE_SWITCH_PARAM,
		RANGE_SWITCH_PARAM,
		COARSE_TUNE_KNOB_PARAM,
		FINE_TUNE_KNOB_PARAM,
		HARM_COUNT_KNOB_PARAM,
		HARM_DENSITY_KNOB_PARAM,
		HARM_STRIDE_KNOB_PARAM,
		HARM_SHIFT_KNOB_PARAM,
		HARM_COUNT_ATTENUVERTER_PARAM,
		HARM_DENSITY_ATTENUVERTER_PARAM,
		LIN_EXP_FM_SWITCH_PARAM,
		HARM_STRIDE_ATTENUVERTER_PARAM,
		HARM_SHIFT_ATTENUVERTER_PARAM,
		SPECTRAL_PIVOT_KNOB_PARAM,
		SPECTRAL_INTENSITY_KNOB_PARAM,
		SPECTRAL_INTENSITY_ATTENUVERTER_PARAM,
		SPECTRAL_TILT_SWITCH_PARAM,
		DRIVE_KNOB_PARAM,
		DRIVE_ATTENUVERTER_PARAM,
		BOOST_FUNDAMENTAL_SWITCH_PARAM,
		OUTPUT_MODE_SWITCH_PARAM,
		PM_ATTENUVERTER_PARAM,
		FM_ATTENUVERTER_PARAM,
		PARAMS_LEN
	};
	enum InputId {
		HARM_COUNT_CV_INPUT,
		HARM_STRIDE_CV_INPUT,
		SPECTRAL_PIVOT_CV_INPUT,
		SPECTRAL_INTENSITY_CV_INPUT,
		DRIVE_CV_INPUT,
		PITCH_INPUT,
		PING_INPUT,
		PM_CV_INPUT,
		FM_CV_INPUT,
		SYNC_INPUT,
		HARM_DENSITY_CV_INPUT,
		HARM_SHIFT_CV_INPUT,
		INPUTS_LEN
	};
	enum OutputId {
		SQUARE_OUTPUT,
		FUNDAMENTAL_OUTPUT,
		ODD_ZERO_DEGREE_OUTPUT,
		EVEN_NINETY_DEGREE_OUTPUT,
		OUTPUTS_LEN
	};
	enum LightId {
		OSCILLATOR_LED_LIGHT_GREEN,
		OSCILLATOR_LED_LIGHT_RED,
		S_LED_1_LIGHT,
		S_LED_2_LIGHT,
		S_LED_3_LIGHT,
		S_LED_4_LIGHT,
		S_LED_5_LIGHT,
		S_LED_6_LIGHT,
		S_LED_7_LIGHT,
		S_LED_8_LIGHT,
		STRIDE_1_LIGHT,
		LIGHTS_LEN
	};

	Loom() {
		config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

		struct CoarseTuneQuantity : ParamQuantity {
			float getDisplayValue() override {
				Loom* module = reinterpret_cast<Loom*>(this->module);
				if (module->lfoMode) {
					// 0.003125Hz (320 second cycle) to 25Hz
					displayMultiplier = LFO_MULTIPLIER;
				} else {
					// 1.25Hz to 10kHz
					displayMultiplier = VCO_MULTIPLIER;
				}
				return ParamQuantity::getDisplayValue();
			}
		};

		struct FineTuneQuantity : ParamQuantity {
			float getDisplayValue() override {
				Loom* module = reinterpret_cast<Loom*>(this->module);
				if (module->lfoMode) {
					// .5x to 2x
					unit = "x";
					displayBase = 2.f;
					displayMultiplier = 1.f;
				} else {
					// -7 to +7 semitones
					unit = " Semitones";
					displayBase = 0.f;
					displayMultiplier = 7.f;
				}
				return ParamQuantity::getDisplayValue();
			}
		};

		struct HarmStrideQuantity : ParamQuantity {
			float getDisplayValue() override {
				Loom* module = reinterpret_cast<Loom*>(this->module);
				auto val = ParamQuantity::getDisplayValue();
				return (module->algo.continuousStrideMode == ContinuousStrideMode::OFF) ? std::round(val) : val;
			}
		};

		// Control Knobs
		configParam<CoarseTuneQuantity>(COARSE_TUNE_KNOB_PARAM, -4.f, 8.965784284662087f, 1.f, "Coarse Tune", " Hz", 2.f);
		configParam<FineTuneQuantity>(FINE_TUNE_KNOB_PARAM, -1.f, 1.f, 0.f, "Fine Tune", " Semitones");
		configParam(HARM_COUNT_KNOB_PARAM, -2.f, 2.f, -2.f, "Harmonic Count", " Partials", 2.f, 16.8, -3.2);
		configParam(HARM_DENSITY_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Density");
		configParam<HarmStrideQuantity>(HARM_STRIDE_KNOB_PARAM, 0.f, 2.f, 0.807354922f, "Harmonic Stride", "x", 2.f, 1.333333333333f, -1.333333333333f);
		configParam(HARM_SHIFT_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Shift");
		configParam(SPECTRAL_PIVOT_KNOB_PARAM, -2.f, 2.f, -2.f, "Spectral Shaping Pivot", "/64", 2.f, 16.8, -3.2);
		configParam(SPECTRAL_INTENSITY_KNOB_PARAM, 0.f, 1.f, 0.f, "Spectral Shaping Intensity");
		configParam(DRIVE_KNOB_PARAM, 0.f, 2.f, 1.f, "Drive");

		// Attenuverters
		configParam(HARM_COUNT_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Harmonic Count CV Attenuverter");
		configParam(HARM_DENSITY_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Harmonic Density CV Attenuverter");
		configParam(HARM_STRIDE_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Harmonic Stride CV Attenuverter");
		configParam(HARM_SHIFT_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Harmonic Shift CV Attenuverter");
		configParam(SPECTRAL_INTENSITY_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Spectral Shaping Intensity CV Attenuverter");
		configParam(DRIVE_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Drive CV Attenuverter");
		configParam(PM_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Phase Modulation CV Attenuverter");
		configParam(FM_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "FM CV Attenuverter");

		// Disable attenuverter randomization
		getParamQuantity(HARM_COUNT_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(HARM_DENSITY_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(HARM_STRIDE_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(HARM_SHIFT_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(SPECTRAL_INTENSITY_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(PM_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(FM_ATTENUVERTER_PARAM)->randomizeEnabled = false;

		// Switches
		configSwitch(CONTINUOUS_STRIDE_SWITCH_PARAM, 0.f, 2.f, 0.f, "Continuous Harmonic Stride", {"Off", "Sync", "Free"});
		configSwitch(RANGE_SWITCH_PARAM, 0.f, 1.f, 1.f, "Oscillator Range", {"LFO", "VCO"});
		configSwitch(LIN_EXP_FM_SWITCH_PARAM, 0.f, 1.f, 0.f, "FM Response", {"Lin", "Exp"});
		configSwitch(SPECTRAL_TILT_SWITCH_PARAM, 0.f, 2.f, 0.f, "Spectral Tilt", {"Lowpass", "Bandpass", "Highpass"});
		configSwitch(BOOST_FUNDAMENTAL_SWITCH_PARAM, 0.f, 1.f, 1.f, "Boost Fundamental", {"Off", "On"});
		configSwitch(OUTPUT_MODE_SWITCH_PARAM, 0.f, 1.f, 0.f, "Output Mode", {"Quadrature", "Odd/Even"});

		// Inputs
		configInput(HARM_COUNT_CV_INPUT, "Harmonic Count CV");
		configInput(HARM_STRIDE_CV_INPUT, "Harmonic Stride CV");
		configInput(SPECTRAL_PIVOT_CV_INPUT, "Spectral Shaping Pivot CV");
		configInput(SPECTRAL_INTENSITY_CV_INPUT, "Spectral Shaping Intensity CV");
		configInput(PITCH_INPUT, "V/Oct Pitch CV");
		configInput(PING_INPUT, "Ping");
		configInput(PM_CV_INPUT, "Phase Modulation CV");
		configInput(FM_CV_INPUT, "FM CV");
		configInput(SYNC_INPUT, "Hard Sync");
		configInput(HARM_DENSITY_CV_INPUT, "Harmonic Density CV");
		configInput(HARM_SHIFT_CV_INPUT, "Harmonic Shift CV");

		// Outputs
		configOutput(SQUARE_OUTPUT, "Square");
		configOutput(FUNDAMENTAL_OUTPUT, "Fundamental (Sine)");
		configOutput(ODD_ZERO_DEGREE_OUTPUT, "Odd / 0°");
		configOutput(EVEN_NINETY_DEGREE_OUTPUT, "Even / 90°");

		lightDivider.setDivision(32);
		pingEnvelope.setRiseFall(650.f, 25.f);
	}

	LoomAlgorithm algo{};
	dsp::ClockDivider lightDivider;
	dsp::ExponentialSlewLimiter pingEnvelope;
	bool lfoMode{false};

	void onReset() override {
		this->algo.oversamplingEnabled = false;
		this->algo.doADAA = true;
		this->algo.doBlep = true;
	}

	json_t* dataToJson() override {
		json_t* rootJ = json_object();
		json_object_set_new(rootJ, "oversample", json_boolean(this->algo.oversamplingEnabled));
		json_object_set_new(rootJ, "doADAA", json_boolean(this->algo.doADAA));
		json_object_set_new(rootJ, "doBlep", json_boolean(this->algo.doBlep));
		return rootJ;
	}

	void dataFromJson(json_t* rootJ) override {
		json_t* oversampleJ = json_object_get(rootJ, "oversample");
		if (oversampleJ) {
			this->algo.oversamplingEnabled = json_boolean_value(oversampleJ);
		}

		json_t* doADAAJ = json_object_get(rootJ, "doADAA");
		if (doADAAJ) {
			this->algo.doADAA = json_boolean_value(doADAAJ);
		}

		json_t* doBlepJ = json_object_get(rootJ, "doBlep");
		if (doBlepJ) {
			this->algo.doBlep = json_boolean_value(doBlepJ);
		}
	}

	void process(const ProcessArgs& args) override {
		// Read and set discrete parameters that aren't interpolated
		this->algo.continuousStrideMode = static_cast<ContinuousStrideMode>(
			(int)params[CONTINUOUS_STRIDE_SWITCH_PARAM].getValue()
		);

		this->algo.boostFund = params[BOOST_FUNDAMENTAL_SWITCH_PARAM].getValue() > .5f;
		this->algo.splitMode = params[OUTPUT_MODE_SWITCH_PARAM].getValue() > .5f;
		this->algo.tilt = (int)params[SPECTRAL_TILT_SWITCH_PARAM].getValue();

		// Get ping envelope so we have it for normalling to unpatched CV inputs
		float pingValue = clamp(inputs[PING_INPUT].getVoltage(), 0.f, 8.f);
		float env = this->pingEnvelope.process(args.sampleTime, pingValue);

		// Read parameters that we'll pack for interpolation, avoiding non-linear calculations
		// such as exponential length or stride scaling
		float length = params[HARM_COUNT_KNOB_PARAM].getValue();
		length += 0.4f * params[HARM_COUNT_ATTENUVERTER_PARAM].getValue() * inputs[HARM_COUNT_CV_INPUT].getNormalVoltage(env);

		float density = params[HARM_DENSITY_KNOB_PARAM].getValue();
		density += 0.2f * params[HARM_DENSITY_ATTENUVERTER_PARAM].getValue() * inputs[HARM_DENSITY_CV_INPUT].getNormalVoltage(env);

		float shift = params[HARM_SHIFT_KNOB_PARAM].getValue();
		shift += 0.2f * params[HARM_SHIFT_ATTENUVERTER_PARAM].getValue() * inputs[HARM_SHIFT_CV_INPUT].getNormalVoltage(env);

		float stride = params[HARM_STRIDE_KNOB_PARAM].getValue();
		stride += 0.2f * params[HARM_STRIDE_ATTENUVERTER_PARAM].getValue() * inputs[HARM_STRIDE_CV_INPUT].getNormalVoltage(env);

		float fmCv = clamp(params[FM_ATTENUVERTER_PARAM].getValue() * inputs[FM_CV_INPUT].getNormalVoltage(env), -5.f, 5.f);
		float coarse = params[COARSE_TUNE_KNOB_PARAM].getValue();
		float fine = params[FINE_TUNE_KNOB_PARAM].getValue();
		float pitchCv = clamp(inputs[PITCH_INPUT].getVoltage(), -2.f, 8.f);

		// Collapse fmCv, coarse, fine, pitchCv into two values, expCv and linCv
		this->lfoMode = params[RANGE_SWITCH_PARAM].getValue() < .5f;
		this->algo.freqMultiplier = this->lfoMode ? LFO_MULTIPLIER : VCO_MULTIPLIER;
		bool expFm = params[LIN_EXP_FM_SWITCH_PARAM].getValue() > .5f;
		fine *= lfoMode ? 1.f : VCO_FINE_TUNE_OCTAVES;
		float expCv = coarse + fine + pitchCv + (expFm ? fmCv : 0.f);
		float linCv = expFm ? 0.f : fmCv * LIN_FM_FACTOR;

		float pivotHarm = params[SPECTRAL_PIVOT_KNOB_PARAM].getValue();
		pivotHarm += 0.4f * inputs[SPECTRAL_PIVOT_CV_INPUT].getVoltage();

		float intensityCv = params[SPECTRAL_INTENSITY_ATTENUVERTER_PARAM].getValue() * .2f * inputs[SPECTRAL_INTENSITY_CV_INPUT].getNormalVoltage(env);
		float intensity = params[SPECTRAL_INTENSITY_KNOB_PARAM].getValue() + intensityCv;

		float driveCv = params[DRIVE_ATTENUVERTER_PARAM].getValue() * .2f * inputs[DRIVE_CV_INPUT].getNormalVoltage(env);
		float drive = params[DRIVE_KNOB_PARAM].getValue() + driveCv;

		float syncValue = inputs[SYNC_INPUT].getVoltage();
		float sinPhaseOffset = params[PM_ATTENUVERTER_PARAM].getValue() * .2f * inputs[PM_CV_INPUT].getNormalVoltage(env);

		// Pack algorithm inputs
		std::array<float_4, 3> algoInputs;
		algoInputs[0] = {length, density, shift, stride};          // Structure section
		algoInputs[1] = {expCv, linCv, sinPhaseOffset, syncValue}; // Pitch/phase influences
		algoInputs[2] = {pivotHarm, intensity, drive, 0.f};        // Shaping section

		// Do stuff
		auto outsPacked = this->algo.process(args, algoInputs)[0];
		outputs[ODD_ZERO_DEGREE_OUTPUT].setVoltage(outsPacked[0]);
		outputs[EVEN_NINETY_DEGREE_OUTPUT].setVoltage(outsPacked[1]);
		outputs[FUNDAMENTAL_OUTPUT].setVoltage(outsPacked[2]);
		outputs[SQUARE_OUTPUT].setVoltage(env);

		if (lightDivider.process()) {
			float lightTime = args.sampleTime * lightDivider.getDivision();
			float oscLight = this->algo.oscLight ? 1.f : -1.f;
			lights[OSCILLATOR_LED_LIGHT_GREEN].setSmoothBrightness(oscLight, lightTime);
			lights[OSCILLATOR_LED_LIGHT_RED].setSmoothBrightness(-oscLight, lightTime);
			lights[STRIDE_1_LIGHT].setSmoothBrightness((float)(this->algo.strideLight), lightTime);
			
			for (int i = 0; i < 8; i++) {
				float brightness = sum_float4(this->algo.ampLights[i]) * .125f;
				lights[S_LED_1_LIGHT + i].setSmoothBrightness(brightness, lightTime);
			}
		}
	}
};

struct LoomWidget : ModuleWidget {
	LoomWidget(Loom* module) {
		setModule(module);
		setPanel(createPanel(asset::plugin(pluginInstance, "res/panels/loom.svg")));

		// Screws
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		// Big knobs
		addParam(createParamCentered<RoundBigBlackKnob>(mm2px(Vec(24.221, 20.552)), module, Loom::COARSE_TUNE_KNOB_PARAM));
		addParam(createParamCentered<RoundBigBlackKnob>(mm2px(Vec(77.059, 22.274)), module, Loom::HARM_COUNT_KNOB_PARAM));

		// Regular knobs
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(10.375, 41.867)), module, Loom::FINE_TUNE_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(81.935, 45.34)), module, Loom::HARM_DENSITY_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(81.935, 62.86)), module, Loom::HARM_SHIFT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(81.935, 81.46)), module, Loom::HARM_STRIDE_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(10.375, 61.879)), module, Loom::SPECTRAL_PIVOT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(10.375, 81.428)), module, Loom::SPECTRAL_INTENSITY_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(39.881, 61.603)), module, Loom::DRIVE_KNOB_PARAM));

		// Trimpots
		addParam(createParamCentered<Trimpot>(mm2px(Vec(56.382, 19.028)), module, Loom::HARM_COUNT_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(67.491, 41.959)), module, Loom::HARM_DENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(67.491, 57.502)), module, Loom::HARM_SHIFT_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(67.491, 73.843)), module, Loom::HARM_STRIDE_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(26.003, 82.006)), module, Loom::SPECTRAL_INTENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(26.003, 61.698)), module, Loom::DRIVE_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(26.073, 41.867)), module, Loom::PM_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(39.881, 41.867)), module, Loom::FM_ATTENUVERTER_PARAM));

		// Vertical switches
		addParam(createParamCentered<CKSS>(mm2px(Vec(8.695, 20.368)), module, Loom::RANGE_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(55.972, 41.279)), module, Loom::BOOST_FUNDAMENTAL_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(55.972, 63.188)), module, Loom::OUTPUT_MODE_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(39.747, 20.552)), module, Loom::LIN_EXP_FM_SWITCH_PARAM));
		addParam(createParamCentered<CKSSThree>(mm2px(Vec(39.777, 81.893)), module, Loom::SPECTRAL_TILT_SWITCH_PARAM));

		// Horizontal switches
		addParam(createParamCentered<CKSSThreeHorizontal>(mm2px(Vec(60.837, 85.522)), module, Loom::CONTINUOUS_STRIDE_SWITCH_PARAM));

		// Inputs
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(29.391, 99.133)), module, Loom::SYNC_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(40.284, 99.133)), module, Loom::HARM_COUNT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(51.177, 99.133)), module, Loom::HARM_STRIDE_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(62.069, 99.133)), module, Loom::SPECTRAL_PIVOT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(62.069, 111.249)), module, Loom::SPECTRAL_INTENSITY_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(29.391, 111.249)), module, Loom::DRIVE_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(18.584, 99.133)), module, Loom::PM_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(18.584, 111.249)), module, Loom::FM_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(40.284, 111.249)), module, Loom::HARM_DENSITY_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(51.177, 111.249)), module, Loom::HARM_SHIFT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(7.585, 99.133)), module, Loom::PITCH_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(7.585, 111.249)), module, Loom::PING_INPUT));

		// Outputs
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(72.976, 99.133)), module, Loom::FUNDAMENTAL_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(83.855, 99.133)), module, Loom::ODD_ZERO_DEGREE_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(72.976, 111.249)), module, Loom::SQUARE_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(83.855, 111.249)), module, Loom::EVEN_NINETY_DEGREE_OUTPUT));

		// Multi-colored LEDs
		addChild(createLightCentered<MediumLight<GreenRedLight>>(mm2px(Vec(33.4, 12.621)), module, Loom::OSCILLATOR_LED_LIGHT_GREEN));

		// Single color LEDs
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(75.410, 74.791)), module, Loom::STRIDE_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(4.233, 73.4)), module, Loom::S_LED_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(8.313, 73.4)), module, Loom::S_LED_2_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(12.393, 73.4)), module, Loom::S_LED_3_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(16.473, 73.4)), module, Loom::S_LED_4_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(20.553, 73.4)), module, Loom::S_LED_5_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(24.633, 73.4)), module, Loom::S_LED_6_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(28.713, 73.4)), module, Loom::S_LED_7_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(32.793, 73.4)), module, Loom::S_LED_8_LIGHT));
	}

	void appendContextMenu(Menu* menu) override {
		Loom* module = dynamic_cast<Loom*>(this->module);

		menu->addChild(new MenuEntry);
		menu->addChild(createMenuLabel("Oversampling"));

		struct OversampleItem : MenuItem {
			Loom* module;
			bool oversample;
			void onAction(const event::Action& e) override {
				module->algo.oversamplingEnabled = oversample;
			}
		};

		{
			OversampleItem* offItem = createMenuItem<OversampleItem>("Off");
			offItem->rightText = CHECKMARK(!(module->algo.oversamplingEnabled));
			offItem->module = module;
			offItem->oversample = false;
			menu->addChild(offItem);

			OversampleItem* onItem = createMenuItem<OversampleItem>("2x");
			onItem->rightText = CHECKMARK(module->algo.oversamplingEnabled);
			onItem->module = module;
			onItem->oversample = true;
			menu->addChild(onItem);
		}

		menu->addChild(new MenuEntry);
		menu->addChild(createMenuLabel("Drive ADAA"));

		struct ADAAItem : MenuItem {
			Loom* module;
			bool doADAA;
			void onAction(const event::Action& e) override {
				module->algo.doADAA = doADAA;
			}
		};

		{
			ADAAItem* offItem = createMenuItem<ADAAItem>("Off");
			offItem->rightText = CHECKMARK(!(module->algo.doADAA));
			offItem->module = module;
			offItem->doADAA = false;
			menu->addChild(offItem);

			ADAAItem* onItem = createMenuItem<ADAAItem>("On");
			onItem->rightText = CHECKMARK(module->algo.doADAA);
			onItem->module = module;
			onItem->doADAA = true;
			menu->addChild(onItem);
		}

		menu->addChild(new MenuEntry);
		menu->addChild(createMenuLabel("MinBLEP Anti-Aliasing"));

		struct BlepItem : MenuItem {
			Loom* module;
			bool doBlep;
			void onAction(const event::Action& e) override {
				module->algo.doBlep = doBlep;
			}
		};

		{
			BlepItem* offItem = createMenuItem<BlepItem>("Off");
			offItem->rightText = CHECKMARK(!(module->algo.doBlep));
			offItem->module = module;
			offItem->doBlep = false;
			menu->addChild(offItem);

			BlepItem* onItem = createMenuItem<BlepItem>("On");
			onItem->rightText = CHECKMARK(module->algo.doBlep);
			onItem->module = module;
			onItem->doBlep = true;
			menu->addChild(onItem);
		}
	}
};


Model* modelLoom = createModel<Loom, LoomWidget>("RigatoniModular-Loom");