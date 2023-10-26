#include "plugin.hpp"
#include <algorithm>
#include <bit>
#include <climits>

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
		HARM_HARM_STRIDE_KNOB_PARAM,
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

		// Control Knobs
		configParam(COARSE_TUNE_KNOB_PARAM, 0.f, 1.f, 0.f, "Coarse Tune");
		configParam(FINE_TUNE_KNOB_PARAM, 0.f, 1.f, 0.f, "Fine Tune");
		configParam(HARM_COUNT_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Count");
		configParam(HARM_DENSITY_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Density");
		configParam(HARM_HARM_STRIDE_KNOB_PARAM, 0.f, 1.f, 0.f, "Harmonic Stride");
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
		this->initialize();
	}

	// All Euclidean pattern bitmaps for each length; inner index is density
	std::array<std::vector<uint64_t>, 64> patternTable;
	std::array<float, 64> harmonicAmplitudes;
	std::array<float, 64> harmonicMultiples;

	// https://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c
	static int uint64_log2(uint64_t n) {
		#define S(k) if (n >= (UINT64_C(1) << k)) { i += k; n >>= k; }
		int i = -(n == 0); S(32); S(16); S(8); S(4); S(2); S(1); return i;
		#undef S
	}

	static uint64_t concatenateBitmaps(uint64_t a, uint64_t b) {
		auto shiftAmount = (b > 1) ? Loom::uint64_log2(b) + 1 : 1;
		return (a << shiftAmount) | b;
	}

	static uint64_t calculateEuclideanBitmap(int length, int density) {
		// Implementation based on https://medium.com/code-music-noise/euclidean-rhythms-391d879494df
		if (length == density) return 0xffffffffffffffff << (64 - length);

		std::vector<uint64_t> ons(density, 1);
		std::vector<uint64_t> offs(length - density, 0);
		while (offs.size() > 1) {
			int numCombinations = std::min(ons.size(), offs.size());
			std::vector<uint64_t> combined{};
			for (int i = 0; i < numCombinations; i++) {
				combined.push_back(Loom::concatenateBitmaps(ons[i], offs[i]));
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
			accum = Loom::concatenateBitmaps(accum, onPattern);
		}
		for (auto &&offPattern : offs) {
			accum = Loom::concatenateBitmaps(accum, offPattern);
		}

		// Make first step of pattern the most significant bit
		return accum << (64 - length);
	}

	void populatePatternTable() {
		for (int length = 1; length <= 64; length++) {
			for (int density = 1; density <= length; density++) {
				patternTable[length - 1].push_back(Loom::calculateEuclideanBitmap(length, density));
			}
		}
	}

	void initialize() {
		harmonicAmplitudes.fill(0.f);
		harmonicMultiples.fill(0.f);
		this->populatePatternTable();
	}

	void process(const ProcessArgs& args) override {
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
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(74.699, 27.816)), module, Loom::HARM_HARM_STRIDE_KNOB_PARAM));
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