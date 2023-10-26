#include "plugin.hpp"
#include <algorithm>
#include <cmath>
#include <climits>

// From VCV Fundamental VCO module
template <typename T>
T sin2pi_pade_05_5_4(T x) {
	x -= 0.5f;
	return (T(-6.283185307) * x + T(33.19863968) * simd::pow(x, 3) - T(32.44191367) * simd::pow(x, 5))
	       / (1 + T(1.296008659) * simd::pow(x, 2) + T(0.7028072946) * simd::pow(x, 4));
}

struct Loom : Module {
	enum ParamId {
		CHARACTER_SWITCH_PARAM,
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
		BOOST_FUNDAMENTAL_SWITCH_PARAM,
		SPECTRAL_PIVOT_KNOB_PARAM,
		SPECTRAL_TILT_KNOB_PARAM,
		SPECTRAL_INTENSITY_KNOB_PARAM,
		HARMONIC_INTENSITY_ATTENUVERTER_PARAM,
		OUTPUT_MODE_SWITCH_PARAM,
		HARMONIC_PIVOT_KNOB_PARAM,
		HARMONIC_TILT_KNOB_PARAM,
		HARMONIC_INTENSITY_KNOB_PARAM,
		SPECTRAL_INTENSITY_ATTENUVERTER_PARAM,
		PM_ATTENUVERTER_PARAM,
		FM_ATTENUVERTER_PARAM,
		PARAMS_LEN
	};
	enum InputId {
		HARM_COUNT_CV_INPUT,
		HARM_STRIDE_CV_INPUT,
		SPECTRAL_PIVOT_CV_INPUT,
		SPECTRAL_TILT_CV_INPUT,
		SPECTRAL_INTENSITY_CV_INPUT,
		PITCH_INPUT,
		PM_CV_INPUT,
		FM_CV_INPUT,
		SYNC_INPUT,
		HARM_DENSITY_CV_INPUT,
		HARM_SHIFT_CV_INPUT,
		WAVESHAPE_PIVOT_CV_INPUT,
		WAVESHAPE_TILT_CV_INPUT,
		WAVESHAPE_INTENSITY_CV_INPUT,
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
		OSCILLATOR_LED_LIGHT,
		S_LED_1_LIGHT,
		S_LED_2_LIGHT,
		S_LED_3_LIGHT,
		S_LED_4_LIGHT,
		S_LED_5_LIGHT,
		H_LED_1_LIGHT,
		H_LED_2_LIGHT,
		H_LED_3_LIGHT,
		H_LED_4_LIGHT,
		H_LED_5_LIGHT,
		STRIDE_1_LIGHT,
		LIGHTS_LEN
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

		// Control Knobs
		configParam<CoarseTuneQuantity>(COARSE_TUNE_KNOB_PARAM, -4.f, 9.f, 1.f, "Coarse Tune", " Hz", 2.f);
		configParam<FineTuneQuantity>(FINE_TUNE_KNOB_PARAM, -1.f, 1.f, 0.f, "Fine Tune", " Semitones");
		configParam(HARM_COUNT_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Count");
		configParam(HARM_DENSITY_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Density");
		configParam(HARM_STRIDE_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Stride");
		configParam(HARM_SHIFT_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Shift");
		configParam(SPECTRAL_PIVOT_KNOB_PARAM, 0.f, 1.f, 0.f, "Spectral Shaping Pivot");
		configParam(SPECTRAL_TILT_KNOB_PARAM, 0.f, 1.f, 0.f, "Spectral Shaping Tilt");
		configParam(SPECTRAL_INTENSITY_KNOB_PARAM, 0.f, 1.f, 0.f, "Spectral Shaping Intensity");
		configParam(HARMONIC_PIVOT_KNOB_PARAM, 0.f, 1.f, 0.f, "Partial Complexity Pivot");
		configParam(HARMONIC_TILT_KNOB_PARAM, 0.f, 1.f, 0.f, "Partial Complexity Tilt");
		configParam(HARMONIC_INTENSITY_KNOB_PARAM, 0.f, 1.f, 0.f, "Partial Complexity Intensity");

		// Attenuverters
		configParam(HARM_COUNT_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Harmonic Count CV Attenuverter");
		configParam(HARM_DENSITY_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Harmonic Density CV Attenuverter");
		configParam(HARM_STRIDE_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Harmonic Stride CV Attenuverter");
		configParam(HARM_SHIFT_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Harmonic Shift CV Attenuverter");
		configParam(HARMONIC_INTENSITY_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Partial Complexity Intensity CV Attenuverter");
		configParam(SPECTRAL_INTENSITY_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Spectral Shaping Intensity CV Attenuverter");
		configParam(PM_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "Phase Modulation CV Attenuverter");
		configParam(FM_ATTENUVERTER_PARAM, -1.f, 1.f, 0.f, "FM CV Attenuverter");

		// Disable attenuverter randomization
		getParamQuantity(HARM_COUNT_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(HARM_DENSITY_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(HARM_STRIDE_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(HARM_SHIFT_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(HARMONIC_INTENSITY_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(SPECTRAL_INTENSITY_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(PM_ATTENUVERTER_PARAM)->randomizeEnabled = false;
		getParamQuantity(FM_ATTENUVERTER_PARAM)->randomizeEnabled = false;

		// Switches
		configSwitch(CHARACTER_SWITCH_PARAM, 0.f, 1.f, 0.f, "Harmonic Distribution Character", {"A", "B"});
		configSwitch(CONTINUOUS_STRIDE_SWITCH_PARAM, 0.f, 1.f, 1.f, "Continuous Harmonic Stride", {"Off", "On"});
		configSwitch(INTERPOLATION_SWITCH_PARAM, 0.f, 1.f, 1.f, "Harmonic Distribution Interpolation", {"Off", "On"});
		configSwitch(RANGE_SWITCH_PARAM, 0.f, 1.f, 1.f, "Oscillator Range", {"LFO", "VCO"});
		configSwitch(LIN_EXP_FM_SWITCH_PARAM, 0.f, 1.f, 0.f, "FM Response", {"Lin", "Exp"});
		configSwitch(BOOST_FUNDAMENTAL_SWITCH_PARAM, 0.f, 1.f, 1.f, "Boost Fundamental", {"Off", "On"});
		configSwitch(OUTPUT_MODE_SWITCH_PARAM, 0.f, 1.f, 1.f, "Output Mode", {"Quadrature", "Odd/Even"});

		// Inputs
		configInput(HARM_COUNT_CV_INPUT, "Harmonic Count CV");
		configInput(HARM_STRIDE_CV_INPUT, "Harmonic Stride CV");
		configInput(SPECTRAL_PIVOT_CV_INPUT, "Spectral Shaping Pivot CV");
		configInput(SPECTRAL_TILT_CV_INPUT, "Spectral Shaping Tilt CV");
		configInput(SPECTRAL_INTENSITY_CV_INPUT, "Spectral Shaping Intensity CV");
		configInput(PITCH_INPUT, "V/Oct Pitch CV");
		configInput(PM_CV_INPUT, "Phase Modulation CV");
		configInput(FM_CV_INPUT, "FM CV");
		configInput(SYNC_INPUT, "Hard Sync");
		configInput(HARM_DENSITY_CV_INPUT, "Harmonic Density CV");
		configInput(HARM_SHIFT_CV_INPUT, "Harmonic Shift CV");
		configInput(WAVESHAPE_PIVOT_CV_INPUT, "Partial Complexity Pivot CV");
		configInput(WAVESHAPE_TILT_CV_INPUT, "Partial Complexity Tilt CV");
		configInput(WAVESHAPE_INTENSITY_CV_INPUT, "Partial Complexity Intensity CV");

		// Outputs
		configOutput(SQUARE_OUTPUT, "Square");
		configOutput(FUNDAMENTAL_OUTPUT, "Fundamental (Sine)");
		configOutput(ODD_ZERO_DEGREE_OUTPUT, "Odd / 0°");
		configOutput(EVEN_NINETY_DEGREE_OUTPUT, "Even / 90°");

		// Set up everything pre-cached/pre-allocated
		this->populatePatternTable();
	}

	// Having this around makes it easier to calculate partial amplitudes from our bitmasks
	static constexpr uint64_t AMP_MASK = 0x8000000000000000;
	static constexpr float LFO_MULTIPLIER = .05f;
	static constexpr float VCO_MULTIPLIER = 20.f;
	static constexpr float LIN_FM_FACTOR = 5.f;
	static constexpr float MAX_FREQ = 10240.f;

	// All Euclidean pattern bitmasks for each length; inner index is density
	std::array<std::vector<uint64_t>, 64> patternTable{};

	// Synthesis parameters
	dsp::MinBlepGenerator<16, 16, float> out1Blep;
	dsp::ClockDivider lightDivider;
	float phaseAccum{0.f};
	bool lfoMode{false};

	// https://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c
	static int uint64_log2(uint64_t n) {
		#define S(k) if (n >= (UINT64_C(1) << k)) { i += k; n >>= k; }
		int i = -(n == 0); S(32); S(16); S(8); S(4); S(2); S(1); return i;
		#undef S
	}

	// Turn patterns into bigger patterns using this one weird trick! CPUs HATE him!
	static uint64_t concatenateBitmasks(uint64_t a, uint64_t b) {
		auto shiftAmount = (b > 1) ? Loom::uint64_log2(b) + 1 : 1;
		return (a << shiftAmount) | b;
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
				patternTable[length - 1].push_back(Loom::calculateEuclideanBitmask(length, density));
			}
		}
	}

	static uint64_t shiftPattern(uint64_t pattern, int shift, int length) {
		return (pattern >> shift) | (pattern << (length - shift));
	}

	int setAmplitudesAndMultiples(
		std::array<float, 64> &amplitudes,
		std::array<float, 64> &multiples,
		float length, // 0-64
		float density, // 0-1
		float stride, // 0-4
		float shift, // 0-1
		bool continuousStride,
		bool interpolate
	) {
		if (!continuousStride) {
			stride = std::round(stride);
		}

		// Simple path
		if (!interpolate) {
			int iLength = (int)std::round(length);
			int iDensity = (int)std::round(density * (iLength - 1));
			int iShift = (int)std::round(shift * iLength) % iLength;

			auto pattern = Loom::shiftPattern(this->patternTable[iLength - 1][iDensity], iShift, iLength);
			for (int i = 0; i < iLength; i++) {
				// TODO: get rid of sawtooth harmonic amplitudes?
				amplitudes[i] = (pattern & AMP_MASK) ? (1.f / (i + 1)) : 0.f;
				multiples[i] = 1.f + i * stride;
				pattern <<= 1;
			}

			return iLength;
		}

		// TODO: interpolation
		return 0;
	}

	static float scaleLength(float normalized) {
		if (normalized >= 0.8f) return 32.f * ((5.f * normalized) - 3.f);
		if (normalized >= 0.6f) return 16.f * ((5.f * normalized) - 2.f);
		if (normalized >= 0.4f) return 8.f * ((5.f * normalized) - 1.f);
		if (normalized >= 0.2f) return 20.f * normalized;
		return 1.f + 15.f * normalized;
	}

	static float scaleStride(float normalized) {
		if (normalized <= 0.666666666667) return 3.f * normalized;
		return (6.f * normalized) - 2.f;
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

	void process(const ProcessArgs& args) override {
		// Read harmonic structure parameters, including attenuverters and CV inputs
		float length = params[HARM_COUNT_KNOB_PARAM].getValue();
		length += 0.2f * params[HARM_COUNT_ATTENUVERTER_PARAM].getValue() * inputs[HARM_COUNT_CV_INPUT].getVoltage();
		length = Loom::scaleLength(clamp(length));

		float density = params[HARM_DENSITY_KNOB_PARAM].getValue();
		density += 0.2f * params[HARM_DENSITY_ATTENUVERTER_PARAM].getValue() * inputs[HARM_DENSITY_CV_INPUT].getVoltage();
		density = clamp(density);

		float stride = params[HARM_STRIDE_KNOB_PARAM].getValue();
		stride += 0.2f * params[HARM_STRIDE_ATTENUVERTER_PARAM].getValue() * inputs[HARM_STRIDE_CV_INPUT].getVoltage();
		stride = Loom::scaleStride(clamp(stride));

		float shift = params[HARM_SHIFT_KNOB_PARAM].getValue();
		shift += 0.2f * params[HARM_SHIFT_ATTENUVERTER_PARAM].getValue() * inputs[HARM_SHIFT_CV_INPUT].getVoltage();
		shift = clamp(shift);

		// Read harmonic structure switches
		bool continuousStride = params[CONTINUOUS_STRIDE_SWITCH_PARAM].getValue() > .5f;
		bool interpolate = params[INTERPOLATION_SWITCH_PARAM].getValue() > .5f;

		std::array<float, 64> harmonicAmplitudes;
		std::array<float, 64> harmonicMultiples;
		int numHarmonics = this->setAmplitudesAndMultiples(
			harmonicAmplitudes, harmonicMultiples, length, density, stride, shift, continuousStride, interpolate
		);

		float out = 0.f;

		// TODO: harmonic complexities array
		// TODO: phase modulation
		// TODO: always update all harmonic phase accumulators every sample
		// TODO: sync all waveforms on changes to continuous stride or VCO/LFO range
		// TODO: sync stride every fundamental cycle if sync stride switch on or continuous stride off

		// Calculate frequency in Hz
		lfoMode = params[RANGE_SWITCH_PARAM].getValue() < .5f;
		bool expFm = params[LIN_EXP_FM_SWITCH_PARAM].getValue() > .5f;
		float fmCv = clamp(params[FM_ATTENUVERTER_PARAM].getValue() * inputs[FM_CV_INPUT].getVoltage(), -5.f, 5.f);
		float coarse = params[COARSE_TUNE_KNOB_PARAM].getValue();
		float fine = params[FINE_TUNE_KNOB_PARAM].getValue();
		float pitchCv = inputs[PITCH_INPUT].getVoltage();
		float freq = Loom::calculateFrequencyHz(coarse, fine, pitchCv, fmCv, expFm, lfoMode);

		this->phaseAccum = this->phaseAccum + freq * args.sampleTime;
		this->phaseAccum -= std::floor(this->phaseAccum);
		for (int i = 0; i < numHarmonics; i++) {
			float harmonicPhaseAccum = this->phaseAccum * harmonicMultiples[i];
			harmonicPhaseAccum -= std::floor(harmonicPhaseAccum);
			out += harmonicAmplitudes[i] * sin2pi_pade_05_5_4(harmonicPhaseAccum);
		}

		// TODO: find a way to normalize output

		out += this->out1Blep.process();
		outputs[ODD_ZERO_DEGREE_OUTPUT].setVoltage(5.f * out);
	}
};

struct CKSSHorizontal : app::SvgSwitch {
	CKSSHorizontal() {
		shadow->opacity = 0.0;
		addFrame(Svg::load(asset::plugin(pluginInstance, "res/components/CKSS_Horizontal_0.svg")));
		addFrame(Svg::load(asset::plugin(pluginInstance, "res/components/CKSS_Horizontal_1.svg")));
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
		addParam(createParamCentered<RoundBigBlackKnob>(mm2px(Vec(21.115, 19.396)), module, Loom::COARSE_TUNE_KNOB_PARAM));

		// Regular knobs
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(21.115, 40.194)), module, Loom::FINE_TUNE_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(40.823, 27.816)), module, Loom::HARM_COUNT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(57.761, 27.816)), module, Loom::HARM_DENSITY_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(74.699, 27.816)), module, Loom::HARM_STRIDE_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(91.637, 27.816)), module, Loom::HARM_SHIFT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(19.853, 61.274)), module, Loom::SPECTRAL_PIVOT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(40.217, 61.274)), module, Loom::SPECTRAL_TILT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(60.581, 61.274)), module, Loom::SPECTRAL_INTENSITY_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(19.853, 81.596)), module, Loom::HARMONIC_PIVOT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(40.217, 81.596)), module, Loom::HARMONIC_TILT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(60.581, 81.627)), module, Loom::HARMONIC_INTENSITY_KNOB_PARAM));

		// Trimpots
		addParam(createParamCentered<Trimpot>(mm2px(Vec(40.823, 44.874)), module, Loom::HARM_COUNT_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(57.761, 44.874)), module, Loom::HARM_DENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(74.699, 44.874)), module, Loom::HARM_STRIDE_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(91.637, 44.874)), module, Loom::HARM_SHIFT_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(78.985, 61.724)), module, Loom::HARMONIC_INTENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(78.985, 81.627)), module, Loom::SPECTRAL_INTENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(5.701, 99.35)), module, Loom::PM_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(16.503, 99.35)), module, Loom::FM_ATTENUVERTER_PARAM));

		// Vertical switches
		addParam(createParamCentered<CKSS>(mm2px(Vec(5.856, 19.212)), module, Loom::RANGE_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(90.872, 61.274)), module, Loom::BOOST_FUNDAMENTAL_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(90.872, 81.596)), module, Loom::OUTPUT_MODE_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(5.856, 40.194)), module, Loom::LIN_EXP_FM_SWITCH_PARAM));

		// Horizontal switches
		addParam(createParamCentered<CKSSHorizontal>(mm2px(Vec(42.9, 15.819)), module, Loom::CHARACTER_SWITCH_PARAM));
		addParam(createParamCentered<CKSSHorizontal>(mm2px(Vec(66.494, 15.819)), module, Loom::CONTINUOUS_STRIDE_SWITCH_PARAM));
		addParam(createParamCentered<CKSSHorizontal>(mm2px(Vec(90.088, 15.819)), module, Loom::INTERPOLATION_SWITCH_PARAM));

		// Inputs
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(27.306, 99.133)), module, Loom::SYNC_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(37.051, 99.133)), module, Loom::HARM_COUNT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(46.795, 99.133)), module, Loom::HARM_STRIDE_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(56.54, 99.133)), module, Loom::SPECTRAL_PIVOT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(66.284, 99.133)), module, Loom::SPECTRAL_TILT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(76.029, 99.133)), module, Loom::SPECTRAL_INTENSITY_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.701, 111.249)), module, Loom::PM_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(16.503, 111.249)), module, Loom::FM_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(27.306, 111.249)), module, Loom::PITCH_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(37.051, 111.249)), module, Loom::HARM_DENSITY_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(46.795, 111.249)), module, Loom::HARM_SHIFT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(56.54, 111.249)), module, Loom::WAVESHAPE_PIVOT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(66.284, 111.249)), module, Loom::WAVESHAPE_TILT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(76.029, 111.249)), module, Loom::WAVESHAPE_INTENSITY_CV_INPUT));

		// Outputs
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(85.773, 99.133)), module, Loom::FUNDAMENTAL_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(95.518, 99.133)), module, Loom::ODD_ZERO_DEGREE_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(85.773, 111.249)), module, Loom::SQUARE_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(95.518, 111.249)), module, Loom::EVEN_NINETY_DEGREE_OUTPUT));

		// Multi-colored LEDs
		addChild(createLightCentered<MediumLight<GreenRedLight>>(mm2px(Vec(11.692, 29.760)), module, Loom::OSCILLATOR_LED_LIGHT));

		// Single color LEDs
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(67.104, 25.882)), module, Loom::STRIDE_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(66.189, 53.278)), module, Loom::S_LED_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(70.247, 53.278)), module, Loom::S_LED_2_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(74.306, 53.278)), module, Loom::S_LED_3_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(78.364, 53.278)), module, Loom::S_LED_4_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(82.422, 53.278)), module, Loom::S_LED_5_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(66.189, 73.305)), module, Loom::H_LED_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(70.247, 73.305)), module, Loom::H_LED_2_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(74.306, 73.305)), module, Loom::H_LED_3_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(78.364, 73.305)), module, Loom::H_LED_4_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(82.422, 73.305)), module, Loom::H_LED_5_LIGHT));
	}
};


Model* modelLoom = createModel<Loom, LoomWidget>("RigatoniModular-Loom");