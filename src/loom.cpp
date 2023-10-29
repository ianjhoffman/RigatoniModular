#include "plugin.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <tuple>
#include <utility>

using simd::float_4;

// From VCV Fundamental VCO module
template <typename T>
T sin2pi_pade_05_5_4(T x) {
	x -= 0.5f;
	return (T(-6.283185307) * x + T(33.19863968) * simd::pow(x, 3) - T(32.44191367) * simd::pow(x, 5))
	       / (1 + T(1.296008659) * simd::pow(x, 2) + T(0.7028072946) * simd::pow(x, 4));
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
		this->phaseAccumulators.fill(0.f);
	}

	// Having this around makes it easier to calculate partial amplitudes from our bitmasks
	static constexpr uint64_t AMP_MASK = 0x8000000000000000;
	static constexpr float LFO_MULTIPLIER = .05f;
	static constexpr float VCO_MULTIPLIER = 20.f;
	static constexpr float LIN_FM_FACTOR = 5.f;
	static constexpr float MAX_FREQ = 10240.f;
	static constexpr float SLOPE_SCALE = 0.0025f;

	// All Euclidean pattern bitmasks for each length; inner index is density, innermost index is step
	// This ends up being a bit wasteful (~200%) but still only uses up 1MB of memory, which I think is fine
	std::array<std::array<std::array<float, 64>, 64>, 64> patternTable{};

	// Synthesis parameters
	std::array<float_4, 16> phaseAccumulators;
	dsp::MinBlepGenerator<16, 16, float> out1Blep, out2Blep, squareBlep, fundBlep;
	dsp::ClockDivider lightDivider;
	float lastFundPhase{0.f}, lastSquare{0.f}, lastSync{0.f};
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

	// Implementation based on https://medium.com/code-music-noise/euclidean-rhythms-391d879494df
	static void calculateEuclideanPattern(std::array<float, 64> &dest, int length, int density) {
		if (length == density) {
			std::fill_n(dest.begin(), length, 1.f);
			return;
		}

		std::vector<uint64_t> ons(density, 1);
		std::vector<uint64_t> offs(length - density, 0);
		while (offs.size() > 1) {
			int numCombinations = std::min(ons.size(), offs.size());
			std::vector<uint64_t> combined{};
			for (int i = 0; i < numCombinations; i++) {
				combined.push_back((ons[i] << Loom::numRelevantBits(offs[i])) | offs[i]);
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

		// Consolidate all patterns, ons then offs, before turning into floats
		ons.insert(ons.end(), offs.begin(), offs.end());
		int idx = 0;
		for (auto pattern : ons) {
			int numSteps = Loom::numRelevantBits(pattern);
			for (int offset = numSteps - 1; offset >= 0; offset--) {
				dest[idx + offset] = simd::ifelse(pattern & 1, 1.f, 0.f);
				pattern >>= 1;
			}
			idx += numSteps;
		}
	}

	void populatePatternTable() {
		for (int length = 1; length <= 64; length++) {
			for (int density = 1; density <= length; density++) {
				Loom::calculateEuclideanPattern(patternTable[length - 1][density - 1], length, density);
			}
		}
	}

	int setAmplitudesAndMultiples(
		std::array<float, 64> &amplitudes,
		std::array<float_4, 16> &multiples,
		float length, // 1-64
		float density, // 0-1
		float stride, // 0-4
		float shift, // 0-1
		bool interpolate
	) {
		// Always calculate all frequency multiples so everything stays aligned
		for (int i = 0; i < 16; i++) {
			int offset = i << 2;
			multiples[i] = {
				1.f + (offset + 0) * stride,
				1.f + (offset + 1) * stride,
				1.f + (offset + 2) * stride,
				1.f + (offset + 3) * stride
			};
		}

		// Simple path
		if (!interpolate) {
			int iLength = (int)std::round(length);
			int iDensity = (int)std::round(density * (iLength - 1));
			int iShift = (int)std::round(shift * (iLength - 1));

			auto pattern = this->patternTable[iLength - 1][iDensity];
			for (int i = 0; i < iLength; i++) {
				amplitudes[i] = pattern[(iShift + i) % iLength];
			}

			return iLength;
		}

		// Interpolation is a lot more complicated but we're essentially finding a point somewhere in the 3D pattern space
		// Length is the "outer" variable in our 2x2x2 fade since it affects the scaling of density and shift
		int iLengthLow = (int)std::floor(length);
		int iLengthHigh = std::min(iLengthLow + 1, 64);
		float lengthFade = length - iLengthLow;

		// Do high length first since it needs to do the extra work of zeroing out the amplitudes
		{
			// Density patterns and fade amount
			float fDensity = density * (iLengthHigh - 1);
			int iDensityLow = (int)std::floor(fDensity);
			int iDensityHigh = iDensityLow + 1;
			float densityFade = fDensity - iDensityLow;

			// Shift values and fade amount
			float fShift = shift * (iLengthHigh - 1);
			int iShiftLow = (int)std::floor(fShift);
			int iShiftHigh = iShiftLow + 1;
			float shiftFade = fShift - iShiftLow;

			float_4 fader = {
				(1.f - densityFade) * (1.f - shiftFade), // low density, low shift
				densityFade * (1.f - shiftFade), // high density, low shift
				(1.f - densityFade) * shiftFade, // low density, high shift
				densityFade * shiftFade // high density, high shift
			};

			// Do blending
			auto lowDens = this->patternTable[iLengthHigh - 1][iDensityLow];
			auto highDens = this->patternTable[iLengthHigh - 1][iDensityHigh];
			for (int i = 0; i < iLengthHigh; i++) {
				int highShiftIdx = (i + iShiftHigh) % iLengthHigh;
				int lowShiftIdx = (i + iShiftLow) % iLengthHigh;
				float_4 highAmp = {
					lowDens[lowShiftIdx], highDens[lowShiftIdx],
					lowDens[highShiftIdx], highDens[highShiftIdx]
				};
				highAmp *= fader * lengthFade;
				amplitudes[i] = highAmp[0] + highAmp[1] + highAmp[2] + highAmp[3];
			}
		}

		// Do low length
		{
			// Density patterns and fade amount
			float fDensity = density * (iLengthLow - 1);
			int iDensityLow = (int)std::floor(fDensity);
			int iDensityHigh = iDensityLow + 1;
			float densityFade = fDensity - iDensityLow;

			// Shift values and fade amount
			float fShift = shift * (iLengthLow - 1);
			int iShiftLow = (int)std::floor(fShift);
			int iShiftHigh = iShiftLow + 1;
			float shiftFade = fShift - iShiftLow;

			float_4 fader = {
				(1.f - densityFade) * (1.f - shiftFade), // low density, low shift
				densityFade * (1.f - shiftFade), // high density, low shift
				(1.f - densityFade) * shiftFade, // low density, high shift
				densityFade * shiftFade // high density, high shift
			};

			// Do blending
			auto lowDens = this->patternTable[iLengthLow - 1][iDensityLow];
			auto highDens = this->patternTable[iLengthLow - 1][iDensityHigh];
			for (int i = 0; i < iLengthLow; i++) {
				int highShiftIdx = (i + iShiftHigh) % iLengthLow;
				int lowShiftIdx = (i + iShiftLow) % iLengthLow;
				float_4 lowAmp = {
					lowDens[lowShiftIdx], highDens[lowShiftIdx],
					lowDens[highShiftIdx], highDens[highShiftIdx]
				};
				lowAmp *= fader * (1.f - lengthFade);
				amplitudes[i] += lowAmp[0] + lowAmp[1] + lowAmp[2] + lowAmp[3];
			}
		}

		return iLengthHigh;
	}

	static std::tuple<float, float, float> calculatePivotSlopes(int tilt) {
		if (tilt == 0) return {0.f, -1.f, 1.f};
		if (tilt == 1) return {-1.5f, -1.5f, 1.5f};
		return {-1.f, 0.f, 1.f};
	}

	static void shapeAmplitudes(
		std::array<float, 64> &amplitudes,
		std::array<float, 5> &ledIndicators,
		float length,
		int tilt, // 0-2 (lowpass, bandpass, highpass)
		float pivot,
		float intensity
	) {
		float cubicPivot = dsp::cubic(pivot);
		float pivotHarm = cubicPivot * (length - 1.f);
		auto slopes = Loom::calculatePivotSlopes(tilt);
		float slopeMultiplier = (intensity < .5f) ? (intensity * 2.f) : (4.f * intensity - 1.f);
		float belowSlope = std::get<0>(slopes) * Loom::SLOPE_SCALE * slopeMultiplier * length;
		float aboveSlope = std::get<1>(slopes) * Loom::SLOPE_SCALE * slopeMultiplier * length;
		float pivotMultiplier = crossfade(1.f, std::get<2>(slopes), clamp(intensity * 2.f));
		for (int i = 0; i < std::ceil(length); i++) {
			float diff = pivotHarm - (float)i;
			float amp = amplitudes[i] * (pivotMultiplier + ((diff > 0) ? belowSlope * diff : aboveSlope * -diff));
			amplitudes[i] = clamp(amp, 0.f, 2.f); // No negative amplitudes
		}

		// LED indicators
		for (int i = 0; i < 5; i++) {
			float diff = pivotHarm - ((length - 1.f) * (0.25f * (float)i));
			ledIndicators[i] = clamp((2.f/3.f) * (pivotMultiplier + ((diff > 0) ? belowSlope * diff : aboveSlope * -diff)));
		}
	}

	static float scaleLength(float normalized) {
		if (normalized >= 0.8f) return 32.f * ((5.f * normalized) - 3.f);
		if (normalized >= 0.6f) return 16.f * ((5.f * normalized) - 2.f);
		if (normalized >= 0.4f) return 8.f * ((5.f * normalized) - 1.f);
		if (normalized >= 0.2f) return 20.f * normalized;
		return 1.f + 15.f * normalized;
	}

	static float scaleStrideKnobValue(float knob) {
		if (knob <= 2.f/3.f) return (3.f/4.f) * knob;
		return (1.5f * knob) - .5f;
	}

	static float calculateFrequencyHz(float coarse, float fine, float pitchCv, float fmCv, bool expFm, bool lfoMode) {
		if (expFm) {
			coarse += fmCv;
		}

		float multiplier = lfoMode ? Loom::LFO_MULTIPLIER : Loom::VCO_MULTIPLIER;
		float fineTuneMultiplier = dsp::exp2_taylor5(fine * (lfoMode ? 1.f : (7.f / 12.f)));
		float freq = multiplier * fineTuneMultiplier * dsp::exp2_taylor5(coarse + pitchCv);

		if (!expFm) {
			freq += fmCv * multiplier * Loom::LIN_FM_FACTOR;
		}

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

		// Actually do partial amplitude/frequency calculations
		std::array<float, 64> harmonicAmplitudes;
		std::array<float_4, 16> harmonicMultiples;
		int numHarmonics = this->setAmplitudesAndMultiples(
			harmonicAmplitudes, harmonicMultiples, length,
			density, stride, shift, this->interpolate
		);

		// Spectral shaping
		float pivot = clamp(params[SPECTRAL_PIVOT_KNOB_PARAM].getValue() + inputs[SPECTRAL_PIVOT_CV_INPUT].getVoltage());
		float tilt = (int)params[SPECTRAL_TILT_SWITCH_PARAM].getValue();
		float intensityCv = params[SPECTRAL_INTENSITY_ATTENUVERTER_PARAM].getValue() * .2f * inputs[SPECTRAL_INTENSITY_CV_INPUT].getVoltage();
		float intensity = clamp(params[SPECTRAL_INTENSITY_KNOB_PARAM].getValue() + intensityCv, -1.f, 1.f);
		std::array<float, 5> shapingLedIndicators;
		Loom::shapeAmplitudes(harmonicAmplitudes, shapingLedIndicators, length, tilt, pivot, intensity);
		
		// Fundamental boosting
		if (params[BOOST_FUNDAMENTAL_SWITCH_PARAM].getValue() > .5f) {
			harmonicAmplitudes[0] = std::max(harmonicAmplitudes[0], .5f);
		}

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

		// Calculate sync
		float syncIn = inputs[SYNC_INPUT].getVoltage();
		float syncCrossingSamplesAgo = -this->lastSync / (syncIn - this->lastSync);
		this->lastSync = syncIn;
		bool sync = (0.f < syncCrossingSamplesAgo) && (syncCrossingSamplesAgo <= 1.f) && (syncIn >= 0.f);

		// Additive synthesis params setup
		bool splitMode = params[OUTPUT_MODE_SWITCH_PARAM].getValue() > .5f;
		float pmCv = clamp(params[PM_ATTENUVERTER_PARAM].getValue() * .2f * inputs[PM_CV_INPUT].getVoltage(), -1.f, 1.f);
		pmCv -= std::floor(pmCv);
		float phaseInc = freq * args.sampleTime;
		float normalSyncPhase = (1.f - syncCrossingSamplesAgo) * phaseInc;
		float fundOut = 0.f, fundOutSync = 0.f, squareOut = 0.f, squareOutSync = 0.f;
		float_4 mainOutsPacked = 0.f;

		// Calculate fundamental sync for non-free continuous stride
		float fundAccum = this->phaseAccumulators[0][0];
		float fundPhaseWrapped = fundAccum + phaseInc - 1.f;
		bool fundSync = fundPhaseWrapped >= 0.f;

		// Some simple built-in band-limiting
		// (not perfect because of drive feature, which can introduce extra harmonics)
		float harmonicMultipleLimit = (args.sampleRate / (2.f * freq)) - 1.f;

		// Calculate this outside the loop to get rid of a bunch of branching
		// 3 types of sync that can introduce discontinuities
		bool doSync = true;
		float syncPhase = fundAccum + phaseInc;
		float anySyncCrossingSamplesAgo = 0.f;
		if (continuousStrideModeChanged || lfoModeChanged) {
			anySyncCrossingSamplesAgo = 0.f;
		} else if (sync) {
			anySyncCrossingSamplesAgo = syncCrossingSamplesAgo;
			syncPhase = normalSyncPhase;
		} else if (continuousStrideMode != ContinuousStrideMode::FREE && fundSync) {
			anySyncCrossingSamplesAgo = fundPhaseWrapped / phaseInc;
		} else {
			doSync = false;
		}

		std::array<float_4, 16> phasesWithoutSync;
		std::array<float_4, 16> phasesWithSync;
		for (int i = 0; i < 16; i++) {
			phasesWithoutSync[i] = this->phaseAccumulators[i] + phaseInc * harmonicMultiples[i];
			phasesWithoutSync[i] -= simd::floor(phasesWithoutSync[i]);
			auto mask = simd::movemaskInverse<float_4>(doSync * 0b1111);
			phasesWithSync[i] = simd::ifelse(mask, harmonicMultiples[i] * syncPhase, phasesWithoutSync[i]);
			phasesWithSync[i] -= simd::floor(phasesWithSync[i]);
			this->phaseAccumulators[i] = phasesWithSync[i];
		}

		for (int i = 0; i < numHarmonics; i++) {
			int f4idx = i >> 2;
			int idx = i & 0b11;
			float harmonicMultiple = harmonicMultiples[f4idx][idx];
			if (harmonicMultiple > harmonicMultipleLimit || harmonicAmplitudes[i] < 0.01f) continue;

			float phaseWithoutSync = phasesWithoutSync[f4idx][idx];
			float phaseWithSync = phasesWithSync[f4idx][idx];

			float_4 packedPhases = {phaseWithoutSync, phaseWithSync, phaseWithoutSync, phaseWithSync};
			float_4 packedPhaseOffsets = {pmCv, pmCv, (.25f + pmCv), (.25f + pmCv)};
			packedPhases += packedPhaseOffsets * harmonicMultiple;
			auto amp = ((float)splitMode + 1.f) * harmonicAmplitudes[i] / (float)(i + 1);
			float_4 splitAmp = simd::ifelse(simd::movemaskInverse<float_4>(splitMode ? 0b1100 : 0), 0.9f, 1.f);

			// [0] = sin without sync, [1] = sin with sync, [2] = cos without sync, [3] = cos with sync
			auto packedVals = sin2pi_pade_05_5_4(packedPhases - simd::floor(packedPhases)) * amp * splitAmp;
			int harmSplitMask = (splitMode && i) ? ((i & 1) ? 0b0011 : 0b1100) : 0b1111;
			mainOutsPacked += simd::ifelse(simd::movemaskInverse<float_4>(harmSplitMask), packedVals, 0.f);
		}

		// TODO: fully implement fundamental and square
		fundAccum = this->phaseAccumulators[0][0];
		fundOut = sin2pi_pade_05_5_4(fundAccum);
		squareOut = simd::ifelse(fundAccum < 0.5f, 1.f, -1.f);

		float driveCv = params[DRIVE_ATTENUVERTER_PARAM].getValue() * .2f * inputs[DRIVE_CV_INPUT].getVoltage();
		float drive = clamp(params[DRIVE_KNOB_PARAM].getValue() + driveCv, 0.f, 2.f);

		// TODO: normalize main outs a bit based on harmonic sum (maybe smart sum based on even vs. odd?)
		mainOutsPacked = 5.f * Loom::drive(mainOutsPacked, drive);

		// BLEP
		if (doSync) {
			float out1Disc = mainOutsPacked[1] - mainOutsPacked[0];
			float out2Disc = mainOutsPacked[3] - mainOutsPacked[2];
			this->out1Blep.insertDiscontinuity(anySyncCrossingSamplesAgo, out1Disc);
			this->out2Blep.insertDiscontinuity(anySyncCrossingSamplesAgo, out2Disc);
		}

		outputs[ODD_ZERO_DEGREE_OUTPUT].setVoltage(mainOutsPacked[1] + this->out1Blep.process());
		outputs[EVEN_NINETY_DEGREE_OUTPUT].setVoltage(mainOutsPacked[3] + this->out2Blep.process());
		outputs[SQUARE_OUTPUT].setVoltage(5.f * (squareOut +  this->squareBlep.process()));
		outputs[FUNDAMENTAL_OUTPUT].setVoltage(5.f * (fundOut + this->fundBlep.process()));

		if (lightDivider.process()) {
			float lightTime = args.sampleTime * lightDivider.getDivision();
			float oscLight = (fundAccum < .5f) ? 1.f : -1.f;
			lights[OSCILLATOR_LED_LIGHT_GREEN].setSmoothBrightness(oscLight, lightTime);
			lights[OSCILLATOR_LED_LIGHT_RED].setSmoothBrightness(-oscLight, lightTime);
			lights[STRIDE_1_LIGHT].setSmoothBrightness(strideIsOne ? 1.f : 0.f, lightTime);
			lights[S_LED_1_LIGHT].setSmoothBrightness(shapingLedIndicators[0], lightTime);
			lights[S_LED_2_LIGHT].setSmoothBrightness(shapingLedIndicators[1], lightTime);
			lights[S_LED_3_LIGHT].setSmoothBrightness(shapingLedIndicators[2], lightTime);
			lights[S_LED_4_LIGHT].setSmoothBrightness(shapingLedIndicators[3], lightTime);
			lights[S_LED_5_LIGHT].setSmoothBrightness(shapingLedIndicators[4], lightTime);
		}

		this->lastSquare = squareOut;
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