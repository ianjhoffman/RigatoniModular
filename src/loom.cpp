#include "plugin.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <tuple>
#include <utility>

using simd::float_4;

// https://web.archive.org/web/20200628195036/http://mooooo.ooo/chebyshev-sine-approximation/
// Thanks to Colin Wallace for doing some great math
const float_4 CHEBYSHEV_COEFFS[6] = {
	float_4(-3.1415926444234477f),   // x
	float_4(2.0261194642649887f),    // x^3
	float_4(-0.5240361513980939f),   // x^5
	float_4(0.0751872634325299f),    // x^7
	float_4(-0.006860187425683514f), // x^9
	float_4(0.000385937753182769f),  // x^11
};

float_4 sin2pi_chebyshev(float_4 x) {
	x = (x * 2.f) - 1.f;
	auto x2 = x * x;
    auto p11 = CHEBYSHEV_COEFFS[5];
    auto p9 = p11 * x2 + CHEBYSHEV_COEFFS[4];
    auto p7 = p9 * x2  + CHEBYSHEV_COEFFS[3];
    auto p5 = p7 * x2  + CHEBYSHEV_COEFFS[2];
    auto p3 = p5 * x2  + CHEBYSHEV_COEFFS[1];
    auto p1 = p3 * x2  + CHEBYSHEV_COEFFS[0];
    return (x - 1.f) * (x + 1.f) * p1 * x;
}

inline float sum_float4(float_4 x) {
	return x[0] + x[1] + x[2] + x[3];
}

struct Loom : Module {
	enum ParamId {
		CONTINUOUS_STRIDE_SWITCH_PARAM,
		INTERPOLATION_SWITCH_PARAM,
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
		STRIDE_1_LIGHT,
		LIGHTS_LEN
	};
	enum ContinuousStrideMode {
		OFF,
		SYNC,
		FREE
	};

	Loom() {
		config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

		struct CoarseTuneQuantity : ParamQuantity {
			float getDisplayValue() override {
				Loom* module = reinterpret_cast<Loom*>(this->module);
				if (module->lfoMode) {
					// 0.003125Hz (320 second cycle) to 25.6Hz
					displayMultiplier = Loom::LFO_MULTIPLIER;
				} else {
					// 1.25Hz to 10.24kHz
					displayMultiplier = Loom::VCO_MULTIPLIER;
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

		struct HarmCountQuantity : ParamQuantity {
			float getDisplayValue() override {
				Loom* module = reinterpret_cast<Loom*>(this->module);
				float out = scaleLength(ParamQuantity::getDisplayValue());
				return module->interpolate ? out : std::round(out);
			}
		};

		struct HarmStrideQuantity : ParamQuantity {
			float getDisplayValue() override {
				Loom* module = reinterpret_cast<Loom*>(this->module);
				float out = 4.f * scaleStrideKnobValue(ParamQuantity::getDisplayValue());
				return module->continuousStride ? out : std::round(out);
			}
		};

		// Control Knobs
		configParam<CoarseTuneQuantity>(COARSE_TUNE_KNOB_PARAM, -4.f, 9.f, 1.f, "Coarse Tune", " Hz", 2.f);
		configParam<FineTuneQuantity>(FINE_TUNE_KNOB_PARAM, -1.f, 1.f, 0.f, "Fine Tune", " Semitones");
		configParam<HarmCountQuantity>(HARM_COUNT_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Count", " Partials");
		configParam(HARM_DENSITY_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Density");
		configParam<HarmStrideQuantity>(HARM_STRIDE_KNOB_PARAM, 0.f, 1.f, 1.f/3.f, "Harmonic Stride", "x");
		configParam(HARM_SHIFT_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Shift");
		configParam(SPECTRAL_PIVOT_KNOB_PARAM, 0.f, 1.f, 0.f, "Spectral Shaping Pivot");
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
		configSwitch(INTERPOLATION_SWITCH_PARAM, 0.f, 1.f, 1.f, "Harmonic Distribution Interpolation", {"Off", "On"});
		configSwitch(RANGE_SWITCH_PARAM, 0.f, 1.f, 1.f, "Oscillator Range", {"LFO", "VCO"});
		configSwitch(LIN_EXP_FM_SWITCH_PARAM, 0.f, 1.f, 0.f, "FM Response", {"Lin", "Exp"});
		configSwitch(SPECTRAL_TILT_SWITCH_PARAM, 0.f, 2.f, 0.f, "Spectral Tilt", {"Lowpass", "Bandpass", "Highpass"});
		configSwitch(BOOST_FUNDAMENTAL_SWITCH_PARAM, 0.f, 1.f, 1.f, "Boost Fundamental", {"Off", "On"});
		configSwitch(OUTPUT_MODE_SWITCH_PARAM, 0.f, 1.f, 1.f, "Output Mode", {"Quadrature", "Odd/Even"});

		// Inputs
		configInput(HARM_COUNT_CV_INPUT, "Harmonic Count CV");
		configInput(HARM_STRIDE_CV_INPUT, "Harmonic Stride CV");
		configInput(SPECTRAL_PIVOT_CV_INPUT, "Spectral Shaping Pivot CV");
		configInput(SPECTRAL_INTENSITY_CV_INPUT, "Spectral Shaping Intensity CV");
		configInput(PITCH_INPUT, "V/Oct Pitch CV");
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

		// Set up everything pre-cached/pre-allocated
		this->populatePatternTable();
		this->populateHarmonicSplitMasks();
		this->calculateBaseAmplitudes();
		this->phaseAccumulators.fill(0.f);
	}

	// Having this around makes it easier to calculate partial amplitudes from our bitmasks
	static constexpr uint64_t AMP_MASK = 0x000000000000000f;
	static constexpr int AMP_SHIFT = 60;
	static constexpr float LFO_MULTIPLIER = .05f;
	static constexpr float VCO_MULTIPLIER = 20.f;
	static constexpr float LIN_FM_FACTOR = 5.f;
	static constexpr float MAX_FREQ = 10240.f;
	static constexpr float SLOPE_SCALE = 0.01f;

	// All Euclidean pattern bitmasks for each length; inner index is density
	std::array<std::vector<uint64_t>, 64> patternTable{};

	// Masks for shifting sequences of each length
	std::array<uint64_t, 64> lengthMasks;

	// Masks for which harmonics go to which output, based on output mode
	std::array<std::array<float_4, 16>, 3> harmonicSplitMasks;

	// Precalculated amplitude multipliers per harmonic to avoid lots of divisions during runtime
	std::array<float_4, 16> baseAmplitudes;

	// Synthesis parameters
	std::array<float_4, 16> phaseAccumulators;
	dsp::MinBlepGenerator<16, 16, float_4> blep;
	dsp::ClockDivider lightDivider;
	float lastSync{0.f};
	ContinuousStrideMode lastContinuousStrideMode{ContinuousStrideMode::OFF};

	// For custom knob displays to read
	bool lfoMode{false};
	bool continuousStride{false};
	bool interpolate{true};

	// https://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c
	static int numRelevantBits(uint64_t n) {
		if (n == 0) return 1;
		#define S(k) if (n >= (UINT64_C(1) << k)) { i += k; n >>= k; }
		int i = -(n == 0); S(32); S(16); S(8); S(4); S(2); S(1); return i + 1;
		#undef S
	}

	inline static uint64_t concatenateBitmasks(uint64_t a, uint64_t b) {
		return (a << Loom::numRelevantBits(b)) | b;
	}

	// Implementation based on https://medium.com/code-music-noise/euclidean-rhythms-391d879494df
	static uint64_t calculateEuclideanBitmask(int length, int density) {
		if (length == density) return 0xffffffffffffffff << (64 - length);

		std::vector<uint64_t> ons(density, 1);
		std::vector<uint64_t> offs(length - density, 0);
		while (offs.size() > 1) {
			int numCombinations = std::min(ons.size(), offs.size());
			std::vector<uint64_t> combined{};
			for (int i = 0; i < numCombinations; i++) {
				combined.push_back(Loom::concatenateBitmasks(ons[i], offs[i]));
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
			accum = Loom::concatenateBitmasks(accum, onPattern);
		}
		for (auto &&offPattern : offs) {
			accum = Loom::concatenateBitmasks(accum, offPattern);
		}

		// Make first step of pattern the most significant bit
		return accum << (64 - length);
	}

	void populatePatternTable() {
		for (int length = 1; length <= 64; length++) {
			for (int density = 1; density <= length; density++) {
				this->patternTable[length - 1].push_back(Loom::calculateEuclideanBitmask(length, density));
			}
			this->lengthMasks[length - 1] = 0xffffffffffffffff << (64 - length);
		}
	}

	inline uint64_t shiftPattern(uint64_t pattern, int shift, int length) {
		auto shifted = ((pattern >> shift) & this->lengthMasks[length - 1]) | (pattern << (length - shift));
		// Swap order of every 4-bit chunk
		uint64_t out = (shifted & 0x8888888888888888) >> 3; // 3 to 0
		out |= (shifted & 0x1111111111111111) << 3; // 0 to 3
		out |= (shifted & 0x4444444444444444) >> 1; // 2 to 1
		out |= (shifted & 0x2222222222222222) << 1; // 1 to 2
		return out;
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

	void setAmplitudesAndMultiples(
		std::array<float_4, 16> &amplitudes,
		std::array<float_4, 16> &multiples,
		float harmonicMultipleLimit,
		float length, // 1-64
		float density, // 0-1
		float stride, // 0-4
		float shift, // 0-1
		bool interpolate
	) {
		// Frequency multiple calculations
		float_4 indices = {0.f, 1.f, 2.f, 3.f};
		for (int i = 0; i < 16; i++) {
			multiples[i] = 1.f + (indices * stride);
			indices += 4.f;
		}

		// Calculate harmonic limit for anti-aliasing, as well as associated mask for shifted patterns
		int harmonicLimit = (harmonicMultipleLimit - 1.f) / std::fmax(stride, 0.1f);
		uint64_t harmonicMask = 0xffffffffffffffff << (64 - clamp(harmonicLimit, 1, 64));

		int iLengthLow = (int)std::floor(length);
		int iLengthHigh = std::min(iLengthLow + 1, 64);
		int numBlocks = (std::min(harmonicLimit, iLengthHigh) + 0b11) >> 2;

		// Simple path
		if (!interpolate) {
			int iLength = (int)std::round(length);
			int lengthIdx = iLength - 1;
			int iDensity = (int)std::round(density * lengthIdx);
			int iShift = (int)std::round(shift * lengthIdx);

			auto pattern = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensity], iShift, iLength);
			for (int i = 0; i < numBlocks; i++) {
				auto ampMask = simd::movemaskInverse<float_4>((pattern >> Loom::AMP_SHIFT) & Loom::AMP_MASK);
				amplitudes[i] = simd::ifelse(ampMask, 1.f, 0.f);
				pattern <<= 4;
			}
		}

		float lengthFade = length - iLengthLow;

		// Do high length first since it needs to do the extra work of zeroing out the amplitudes
		{
			int lengthIdx = iLengthHigh - 1;

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
			auto lowDensLowShift = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensityLow], iShiftLow, iLengthHigh);
			auto lowDensHighShift = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensityLow], iShiftHigh, iLengthHigh);
			auto highDensLowShift = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensityHigh], iShiftLow, iLengthHigh);
			auto highDensHighShift = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensityHigh], iShiftHigh, iLengthHigh);
			auto shiftAmt = Loom::AMP_SHIFT;
			for (int i = 0; i < numBlocks; i++) {
				auto lowLow = simd::movemaskInverse<float_4>((lowDensLowShift >> shiftAmt) & Loom::AMP_MASK);
				auto lowHigh = simd::movemaskInverse<float_4>((lowDensHighShift >> shiftAmt) & Loom::AMP_MASK);
				auto highLow = simd::movemaskInverse<float_4>((highDensLowShift >> shiftAmt) & Loom::AMP_MASK);
				auto highHigh = simd::movemaskInverse<float_4>((highDensHighShift >> shiftAmt) & Loom::AMP_MASK);

				// Wish I could do a matrix multiplication here
				amplitudes[i] += simd::ifelse(lowLow, fader[0], 0.f);
				amplitudes[i] += simd::ifelse(lowHigh, fader[1], 0.f);
				amplitudes[i] += simd::ifelse(highLow, fader[2], 0.f);
				amplitudes[i] += simd::ifelse(highHigh, fader[3], 0.f);

				shiftAmt -= 4;
			}
		}

		// Do low length
		{
			int lengthIdx = iLengthLow - 1;

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
			fader *= (1.f - lengthFade);

			// Do blending
			auto lowDensLowShift = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensityLow], iShiftLow, iLengthLow);
			auto lowDensHighShift = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensityLow], iShiftHigh, iLengthLow);
			auto highDensLowShift = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensityHigh], iShiftLow, iLengthLow);
			auto highDensHighShift = harmonicMask & this->shiftPattern(this->patternTable[lengthIdx][iDensityHigh], iShiftHigh, iLengthLow);
			auto shiftAmt = Loom::AMP_SHIFT;
			for (int i = 0; i < numBlocks; i++) {
				auto lowLow = simd::movemaskInverse<float_4>((lowDensLowShift >> shiftAmt) & Loom::AMP_MASK);
				auto lowHigh = simd::movemaskInverse<float_4>((lowDensHighShift >> shiftAmt) & Loom::AMP_MASK);
				auto highLow = simd::movemaskInverse<float_4>((highDensLowShift >> shiftAmt) & Loom::AMP_MASK);
				auto highHigh = simd::movemaskInverse<float_4>((highDensHighShift >> shiftAmt) & Loom::AMP_MASK);

				// Wish I could do a matrix multiplication here
				amplitudes[i] += simd::ifelse(lowLow, fader[0], 0.f);
				amplitudes[i] += simd::ifelse(lowHigh, fader[1], 0.f);
				amplitudes[i] += simd::ifelse(highLow, fader[2], 0.f);
				amplitudes[i] += simd::ifelse(highHigh, fader[3], 0.f);

				shiftAmt -= 4;
			}
		}
	}

	inline static std::pair<float, float> calculatePivotSlopes(int tilt) {
		return (tilt == 0) ? std::make_pair(0.1f, -2.f)
			: ((tilt == 1) ? std::make_pair(-3.f, -3.f)
			: std::make_pair(-2.f, 0.1f));
	}

	static void shapeAmplitudes(
		std::array<float_4, 16> &amplitudes,
		std::array<float, 5> &ledIndicators,
		float length,
		int tilt, // 0-2 (lowpass, bandpass, highpass)
		float pivot,
		float intensity
	) {
		float cubicPivot = dsp::cubic(pivot);
		float pivotHarm = cubicPivot * 63;
		auto slopes = Loom::calculatePivotSlopes(tilt);
		float slopeMultiplier = (intensity < .5f) ? (intensity * 8.f) : (16.f * intensity - 4.f);
		float belowSlope = std::get<0>(slopes) * Loom::SLOPE_SCALE * slopeMultiplier;
		float aboveSlope = std::get<1>(slopes) * Loom::SLOPE_SCALE * slopeMultiplier;
		float pivotBase = crossfade(0.f, (tilt == 1) ? .15f : 0.f, clamp(intensity * 2.f));
		float_4 indices = {0.f, 1.f, 2.f, 3.f};
		for (int i = 0; i < 16; i++) {
			float_4 diffs = pivotHarm - indices;
			float_4 slopes = simd::ifelse(diffs > 0, belowSlope, aboveSlope);
			auto toAdd = pivotBase + slopes * simd::abs(diffs);
			amplitudes[i] = clamp(amplitudes[i] + toAdd, 0.f, 1.5f); // No negative amplitudes
			indices += 4.f;
		}

		// LED indicators
		for (int i = 0; i < 5; i++) {
			// Always use max harmonic length for LED display
			float diff = pivotHarm - (63 * (0.25f * (float)i));
			ledIndicators[i] = clamp(4.f * (pivotBase + ((diff > 0) ? belowSlope * diff : aboveSlope * -diff)));
		}
	}

	static float scaleLength(float normalized) {
		auto pow = (normalized <= 0.2f) ? (10.f * normalized) : (5.f * normalized + 1.f);
		return clamp(dsp::exp2_taylor5(pow), 1.f, 64.f);
	}

	static float scaleStrideKnobValue(float knob) {
		return (knob <= 2.f/3.f) ? ((3.f/4.f) * knob) : ((1.5f * knob) - .5f);
	}

	static float calculateFrequencyHz(float coarse, float fine, float pitchCv, float fmCv, bool expFm, bool lfoMode) {
		coarse += expFm ? fmCv : 0.f;
		float multiplier = lfoMode ? Loom::LFO_MULTIPLIER : Loom::VCO_MULTIPLIER;
		float fineTuneMultiplier = dsp::exp2_taylor5(fine * (lfoMode ? 1.f : (7.f / 12.f)));
		float freq = multiplier * fineTuneMultiplier * dsp::exp2_taylor5(coarse + pitchCv);
		freq += expFm ? 0.f : fmCv * multiplier * Loom::LIN_FM_FACTOR;

		// Clamp to max frequency, allow negative frequency for thru-zero linear FM
		return std::min(freq, Loom::MAX_FREQ);
	}

	// Soft clipping (quadratic) drive with a small amount of wave folding at the extreme
	static float_4 drive(float_4 in, float drive) {
		// Don't hit the folder as hard, up to 12 o'clock on the drive knob should almost never clip
		in *= drive * .5f;
		auto overThresh = simd::abs(in) - .9f;
		auto underThreshMask = overThresh < 0.f;
		auto subAmt = simd::clamp(overThresh * overThresh * drive, 0.f, overThresh * 1.25f);
		auto driven = simd::ifelse(in < 0.f, in + subAmt, in - subAmt);
		return simd::ifelse(underThreshMask, in, driven);
	}

	void process(const ProcessArgs& args) override {
		// Read harmonic structure switches
		float continuousStrideParam = params[CONTINUOUS_STRIDE_SWITCH_PARAM].getValue();
		auto continuousStrideMode = (continuousStrideParam < .5f) ? ContinuousStrideMode::OFF :
			((continuousStrideParam < 1.5f) ? ContinuousStrideMode::SYNC : ContinuousStrideMode::FREE);
		bool continuousStrideModeChanged = continuousStrideMode != this->lastContinuousStrideMode;
		this->lastContinuousStrideMode = continuousStrideMode;
		this->continuousStride = continuousStrideMode != ContinuousStrideMode::OFF;
		this->interpolate = params[INTERPOLATION_SWITCH_PARAM].getValue() > .5f;

		// Read harmonic structure parameters, including attenuverters and CV inputs
		float length = params[HARM_COUNT_KNOB_PARAM].getValue();
		length += 0.2f * params[HARM_COUNT_ATTENUVERTER_PARAM].getValue() * inputs[HARM_COUNT_CV_INPUT].getVoltage();
		length = Loom::scaleLength(clamp(length));

		float density = params[HARM_DENSITY_KNOB_PARAM].getValue();
		density += 0.2f * params[HARM_DENSITY_ATTENUVERTER_PARAM].getValue() * inputs[HARM_DENSITY_CV_INPUT].getVoltage();
		density = clamp(density);

		float strideScaled = Loom::scaleStrideKnobValue(params[HARM_STRIDE_KNOB_PARAM].getValue());
		float stride = 4.f * clamp(
			strideScaled + 0.1f * params[HARM_STRIDE_ATTENUVERTER_PARAM].getValue() * inputs[HARM_STRIDE_CV_INPUT].getVoltage()
		);
		if (!this->continuousStride) {
			stride = std::round(stride);
		}
		bool strideIsOne = std::abs(stride - 1.f) < .01f; // Reasonable epsilon

		float shift = params[HARM_SHIFT_KNOB_PARAM].getValue();
		shift += 0.2f * params[HARM_SHIFT_ATTENUVERTER_PARAM].getValue() * inputs[HARM_SHIFT_CV_INPUT].getVoltage();
		shift = clamp(shift);

		// Calculate fundamental frequency in Hz
		bool newLfoMode = params[RANGE_SWITCH_PARAM].getValue() < .5f;
		bool lfoModeChanged = newLfoMode != this->lfoMode;
		this->lfoMode = newLfoMode;
		bool expFm = params[LIN_EXP_FM_SWITCH_PARAM].getValue() > .5f;
		float fmCv = clamp(params[FM_ATTENUVERTER_PARAM].getValue() * inputs[FM_CV_INPUT].getVoltage(), -5.f, 5.f);
		float coarse = params[COARSE_TUNE_KNOB_PARAM].getValue();
		float fine = params[FINE_TUNE_KNOB_PARAM].getValue();
		float pitchCv = inputs[PITCH_INPUT].getVoltage();
		float freq = Loom::calculateFrequencyHz(coarse, fine, pitchCv, fmCv, expFm, this->lfoMode);

		// Some simple built-in band-limiting (not perfect because of
		// drive, fm, pm, all of which can introduce extra harmonics)
		auto freq2Recip = simd::rcp(2.f * freq);
		float harmonicMultipleLimit = args.sampleRate * freq2Recip[0] - 1.f;

		// Actually do partial amplitude/frequency calculations
		std::array<float_4, 16> harmonicAmplitudes{};
		std::array<float_4, 16> harmonicMultiples{};
		this->setAmplitudesAndMultiples(
			harmonicAmplitudes, harmonicMultiples, harmonicMultipleLimit, length, density, stride, shift, this->interpolate
		);

		// Spectral shaping
		float pivot = clamp(params[SPECTRAL_PIVOT_KNOB_PARAM].getValue() + inputs[SPECTRAL_PIVOT_CV_INPUT].getVoltage());
		float tilt = (int)params[SPECTRAL_TILT_SWITCH_PARAM].getValue();
		float intensityCv = params[SPECTRAL_INTENSITY_ATTENUVERTER_PARAM].getValue() * .2f * inputs[SPECTRAL_INTENSITY_CV_INPUT].getVoltage();
		float intensity = clamp(params[SPECTRAL_INTENSITY_KNOB_PARAM].getValue() + intensityCv, -1.f, 1.f);
		std::array<float, 5> shapingLedIndicators;
		Loom::shapeAmplitudes(harmonicAmplitudes, shapingLedIndicators, length, tilt, pivot, intensity);
		
		// Fundamental boosting
		bool boostFund = params[BOOST_FUNDAMENTAL_SWITCH_PARAM].getValue() > .5f;
		harmonicAmplitudes[0][0] = std::max(harmonicAmplitudes[0][0], boostFund ? .5f : 0.f);

		// Calculate sync
		float syncIn = inputs[SYNC_INPUT].getVoltage();
		float syncCrossing = -this->lastSync / (syncIn - this->lastSync);
		this->lastSync = syncIn;
		bool sync = (0.f < syncCrossing) && (syncCrossing <= 1.f) && (syncIn >= 0.f);

		// Additive synthesis params setup
		bool splitMode = params[OUTPUT_MODE_SWITCH_PARAM].getValue() > .5f;
		float phaseOffset = clamp(params[PM_ATTENUVERTER_PARAM].getValue() * .2f * inputs[PM_CV_INPUT].getVoltage(), -1.f, 1.f);
		phaseOffset -= std::floor(phaseOffset);
		float cosPhaseOffset = phaseOffset + 0.25f;
		float phaseInc = freq * args.sampleTime;
		float normalSyncPhase = (1.f - syncCrossing) * phaseInc;

		// Calculate fundamental sync for non-free continuous stride
		float fundAccum = this->phaseAccumulators[0][0];
		float fundPhaseWrapped = fundAccum + phaseInc - 1.f;
		bool fundSync = fundPhaseWrapped >= 0.f;

		// Calculate this outside the loop to get rid of a bunch of branching
		// 3 types of sync that can introduce discontinuities
		bool doSync = true;
		float syncPhase = fundAccum + phaseInc;
		float minBlepP = 0.f;
		if (continuousStrideModeChanged || lfoModeChanged) {
			minBlepP = 0.f;
		} else if (sync) {
			minBlepP = syncCrossing - 1.f;
			syncPhase = normalSyncPhase;
		} else if (continuousStrideMode != ContinuousStrideMode::FREE && fundSync) {
			minBlepP = (fundPhaseWrapped / phaseInc) - 1.f;
		} else {
			doSync = false;
		}

		std::array<float_4, 16> out1WithoutSync;
		std::array<float_4, 16> out1WithSync;
		std::array<float_4, 16> out2WithoutSync;
		std::array<float_4, 16> out2WithSync;
		int oddHarmSplitMaskIdx = splitMode ? 1 : 0;
		int evenHarmSplitMaskIdx = splitMode ? 2 : 0;
		auto syncMask = simd::movemaskInverse<float_4>(doSync ? 0b1111 : 0b0000);
		auto cosPhaseMask = simd::movemaskInverse<float_4>(splitMode ? 0b0000 : 0b1111);
		float_4 out1WithoutSyncSum{}, out1WithSyncSum{}, out2WithoutSyncSum{}, out2WithSyncSum{};
		float_4 ampMult = (float)splitMode + 1.f;
		for (int i = 0; i < 16; i++) {
			float_4 overallAmplitude = harmonicAmplitudes[i] * this->baseAmplitudes[i] * ampMult;

			out1WithoutSync[i] = this->phaseAccumulators[i] + (phaseInc * harmonicMultiples[i]);
			out1WithoutSync[i] -= simd::floor(out1WithoutSync[i]);
			out1WithSync[i] = simd::ifelse(syncMask, harmonicMultiples[i] * syncPhase, out1WithoutSync[i]);
			out1WithSync[i] -= simd::floor(out1WithSync[i]);
			this->phaseAccumulators[i] = out1WithSync[i]; // store before doing PM

			// Phase modulation
			out2WithoutSync[i] = out1WithoutSync[i] + (cosPhaseOffset * harmonicMultiples[i]);
			out2WithoutSync[i] -= simd::floor(out2WithoutSync[i]);

			out2WithSync[i] = out1WithSync[i] + (cosPhaseOffset * harmonicMultiples[i]);
			out2WithSync[i] -= simd::floor(out2WithSync[i]);

			out1WithoutSync[i] += (phaseOffset * harmonicMultiples[i]);
			out1WithoutSync[i] -= simd::floor(out1WithoutSync[i]);

			out1WithSync[i] += (phaseOffset * harmonicMultiples[i]);
			out1WithSync[i] -= simd::floor(out1WithSync[i]);

			// Output 2 doesn't use cosine phase partials in odd/even split mode
			out2WithoutSync[i] = simd::ifelse(cosPhaseMask, out2WithoutSync[i], out1WithoutSync[i]);
			out2WithSync[i] = simd::ifelse(cosPhaseMask, out2WithSync[i], out1WithSync[i]);

			// Calculate sin/cos with amplitudes, sum with output
			out1WithoutSyncSum += simd::ifelse(
				this->harmonicSplitMasks[oddHarmSplitMaskIdx][i],
				sin2pi_chebyshev(out1WithoutSync[i]) * overallAmplitude,
				0.f
			);

			out1WithSyncSum += simd::ifelse(
				this->harmonicSplitMasks[oddHarmSplitMaskIdx][i],
				sin2pi_chebyshev(out1WithSync[i]) * overallAmplitude,
				0.f
			);

			out2WithoutSyncSum += simd::ifelse(
				this->harmonicSplitMasks[evenHarmSplitMaskIdx][i],
				sin2pi_chebyshev(out2WithoutSync[i]) * overallAmplitude,
				0.f
			);

			out2WithSyncSum += simd::ifelse(
				this->harmonicSplitMasks[evenHarmSplitMaskIdx][i],
				sin2pi_chebyshev(out2WithSync[i]) * overallAmplitude,
				0.f
			);
		}

		// Collapse the remaining 4 harmonic accumulators for each output
		float_4 mainOutsPacked = {out1WithoutSyncSum[0], out1WithSyncSum[0], out2WithoutSyncSum[0], out2WithSyncSum[0]};
		mainOutsPacked += {out1WithoutSyncSum[1], out1WithSyncSum[1], out2WithoutSyncSum[1], out2WithSyncSum[1]};
		mainOutsPacked += {out1WithoutSyncSum[2], out1WithSyncSum[2], out2WithoutSyncSum[2], out2WithSyncSum[2]};
		mainOutsPacked += {out1WithoutSyncSum[3], out1WithSyncSum[3], out2WithoutSyncSum[3], out2WithSyncSum[3]};

		// TODO: fully implement fundamental and square
		auto fundPhaseWithoutSync = out1WithoutSync[0][0];
		auto fundPhaseWithSync = out1WithSync[0][0];
		float fundOutWithoutSync = simd::sin(2.f * M_PI * fundPhaseWithoutSync);
		float fundOutWithSync = simd::sin(2.f * M_PI * fundPhaseWithSync);
		auto squareOutWithoutSync = simd::ifelse(fundPhaseWithoutSync < 0.5f, 1.f, -1.f);
		auto squareOutWithSync = simd::ifelse(fundPhaseWithSync < 0.5f, 1.f, -1.f);

		float driveCv = params[DRIVE_ATTENUVERTER_PARAM].getValue() * .2f * inputs[DRIVE_CV_INPUT].getVoltage();
		float drive = clamp(params[DRIVE_KNOB_PARAM].getValue() + driveCv, 0.f, 2.f);
		mainOutsPacked = Loom::drive(mainOutsPacked, drive);

		// BLEP
		if (doSync) {
			float_4 discontinuities = {
				mainOutsPacked[1] - mainOutsPacked[0],
				mainOutsPacked[3] - mainOutsPacked[2],
				fundOutWithSync - fundOutWithoutSync,
				squareOutWithSync - squareOutWithoutSync
			};
			this->blep.insertDiscontinuity(minBlepP, discontinuities);
		}

		float_4 outsPacked  = {mainOutsPacked[1], mainOutsPacked[3], fundOutWithSync, squareOutWithSync};
		outsPacked = 5.f * (outsPacked + this->blep.process());

		outputs[ODD_ZERO_DEGREE_OUTPUT].setVoltage(outsPacked[0]);
		outputs[EVEN_NINETY_DEGREE_OUTPUT].setVoltage(outsPacked[1]);
		outputs[FUNDAMENTAL_OUTPUT].setVoltage(outsPacked[2]);
		outputs[SQUARE_OUTPUT].setVoltage(outsPacked[3]);

		if (lightDivider.process()) {
			float lightTime = args.sampleTime * lightDivider.getDivision();
			float oscLight = (fundAccum < .5f) ? 1.f : -1.f;
			lights[OSCILLATOR_LED_LIGHT_GREEN].setSmoothBrightness(oscLight, lightTime);
			lights[OSCILLATOR_LED_LIGHT_RED].setSmoothBrightness(-oscLight, lightTime);
			lights[STRIDE_1_LIGHT].setSmoothBrightness(strideIsOne ? 1.f : 0.f, lightTime);
			// TODO - 8 LEDs, each showing 8 partials' amplitudes
			lights[S_LED_1_LIGHT].setSmoothBrightness(shapingLedIndicators[0], lightTime);
			lights[S_LED_2_LIGHT].setSmoothBrightness(shapingLedIndicators[1], lightTime);
			lights[S_LED_3_LIGHT].setSmoothBrightness(shapingLedIndicators[2], lightTime);
			lights[S_LED_4_LIGHT].setSmoothBrightness(shapingLedIndicators[3], lightTime);
			lights[S_LED_5_LIGHT].setSmoothBrightness(shapingLedIndicators[4], lightTime);
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
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(10.375, 37.351)), module, Loom::FINE_TUNE_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(81.292, 45.34)), module, Loom::HARM_DENSITY_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(81.292, 62.86)), module, Loom::HARM_SHIFT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(81.292, 81.46)), module, Loom::HARM_STRIDE_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(10.375, 59.763)), module, Loom::SPECTRAL_PIVOT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(10.375, 81.428)), module, Loom::SPECTRAL_INTENSITY_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(39.777, 68.589)), module, Loom::DRIVE_KNOB_PARAM));

		// Trimpots
		addParam(createParamCentered<Trimpot>(mm2px(Vec(58.498, 15.157)), module, Loom::HARM_COUNT_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(66.848, 41.959)), module, Loom::HARM_DENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(66.848, 56.443)), module, Loom::HARM_SHIFT_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(66.848, 71.726)), module, Loom::HARM_STRIDE_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(26.003, 82.006)), module, Loom::SPECTRAL_INTENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(39.773, 86.016)), module, Loom::DRIVE_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(26.215, 41.127)), module, Loom::PM_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(39.861, 37.845)), module, Loom::FM_ATTENUVERTER_PARAM));

		// Vertical switches
		addParam(createParamCentered<CKSS>(mm2px(Vec(5.787, 20.368)), module, Loom::RANGE_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(55.443, 49.195)), module, Loom::BOOST_FUNDAMENTAL_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(55.443, 68.458)), module, Loom::OUTPUT_MODE_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(42.119, 20.552)), module, Loom::LIN_EXP_FM_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(55.443, 30.088)), module, Loom::INTERPOLATION_SWITCH_PARAM));

		// Horizontal switches
		addParam(createParamCentered<CKSSThreeHorizontal>(mm2px(Vec(60.837, 85.522)), module, Loom::CONTINUOUS_STRIDE_SWITCH_PARAM));
		addParam(createParamCentered<CKSSThreeHorizontal>(mm2px(Vec(31.145, 55.811)), module, Loom::SPECTRAL_TILT_SWITCH_PARAM));

		// Inputs
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(29.580, 99.133)), module, Loom::SYNC_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(40.473, 99.133)), module, Loom::HARM_COUNT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(51.366, 99.133)), module, Loom::HARM_STRIDE_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(62.258, 99.133)), module, Loom::SPECTRAL_PIVOT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(62.258, 111.249)), module, Loom::SPECTRAL_INTENSITY_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(29.580, 111.249)), module, Loom::DRIVE_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(18.773, 99.133)), module, Loom::PM_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(18.773, 111.249)), module, Loom::FM_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(40.473, 111.249)), module, Loom::HARM_DENSITY_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(51.366, 111.249)), module, Loom::HARM_SHIFT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(7.774, 105.553)), module, Loom::PITCH_INPUT));

		// Outputs
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(73.165, 99.133)), module, Loom::FUNDAMENTAL_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(84.044, 99.133)), module, Loom::ODD_ZERO_DEGREE_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(73.165, 111.249)), module, Loom::SQUARE_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(84.044, 111.249)), module, Loom::EVEN_NINETY_DEGREE_OUTPUT));

		// Multi-colored LEDs
		addChild(createLightCentered<MediumLight<GreenRedLight>>(mm2px(Vec(13.583, 13.071)), module, Loom::OSCILLATOR_LED_LIGHT_GREEN));

		// Single color LEDs
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(71.907, 78.750)), module, Loom::STRIDE_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(12.515, 72.45)), module, Loom::S_LED_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(16.574, 72.45)), module, Loom::S_LED_2_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(20.632, 72.45)), module, Loom::S_LED_3_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(24.690, 72.45)), module, Loom::S_LED_4_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(28.749, 72.45)), module, Loom::S_LED_5_LIGHT));
	}
};


Model* modelLoom = createModel<Loom, LoomWidget>("RigatoniModular-Loom");