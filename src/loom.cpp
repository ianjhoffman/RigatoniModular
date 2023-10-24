#include "plugin.hpp"


struct Loom : Module {
	enum ParamId {
		CHARACTER_SWITCH_PARAM,
		CONTINUOUS_STRIDE_SWITCH_PARAM,
		INTERPOLATION_SWITCH_PARAM,
		RANGE_SWITCH_PARAM,
		COARSE_TUNE_KNOB_PARAM,
		HARM_COUNT_KNOB_PARAM,
		HARM_DENSITY_KNOB_PARAM,
		HARM_HARM_STRIDE_KNOB_PARAM,
		HARM_SHIFT_KNOB_PARAM,
		HARM_COUNT_ATTENUVERTER_PARAM,
		HARM_DENSITY_ATTENUVERTER_PARAM,
		LIN_EXP_FM_SWITCH_PARAM,
		HARM_STRIDE_ATTENUVERTER_PARAM,
		HARM_SHIFT_ATTENUVERTER_PARAM,
		HARMONIC_CURVE_SWITCH_PARAM,
		SPECTRAL_PIVOT_KNOB_PARAM,
		SPECTRAL_TILT_KNOB_PARAM,
		SPECTRAL_INTENSITY_KNOB_PARAM,
		HARMONIC_INTENSITY_ATTENUVERTER_PARAM,
		SPECTRAL_CURVE_SWITCH_PARAM,
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
		SQUARE_OUT_OUTPUT,
		MAIN_OUT_OUTPUT,
		SUB_OUT_OUTPUT,
		QUADRATURE_OUT_OUTPUT,
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

		// 2-Value Switches
		configSwitch(CHARACTER_SWITCH_PARAM, 0.f, 1.f, 0.f, "Harmonic Distribution Character", {"A", "B"});
		configSwitch(CONTINUOUS_STRIDE_SWITCH_PARAM, 0.f, 1.f, 1.f, "Continuous Harmonic Stride", {"Off", "On"});
		configSwitch(INTERPOLATION_SWITCH_PARAM, 0.f, 1.f, 1.f, "Harmonic Distribution Interpolation", {"Off", "On"});
		configSwitch(RANGE_SWITCH_PARAM, 0.f, 1.f, 1.f, "Oscillator Range", {"LFO", "VCO"});
		configSwitch(LIN_EXP_FM_SWITCH_PARAM, 0.f, 1.f, 0.f, "FM Response", {"Lin", "Exp"});

		// 3-Value Switches
		configSwitch(HARMONIC_CURVE_SWITCH_PARAM, 0.f, 2.f, 1.f, "Partial Complexity Curve", {"Log", "Lin", "Exp"});
		configSwitch(SPECTRAL_CURVE_SWITCH_PARAM, 0.f, 2.f, 1.f, "Spectral Shaping Curve", {"Log", "Lin", "Exp"});

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
		configOutput(SQUARE_OUT_OUTPUT, "Square");
		configOutput(MAIN_OUT_OUTPUT, "Main");
		configOutput(SUB_OUT_OUTPUT, "Sub");
		configOutput(QUADRATURE_OUT_OUTPUT, "Quadrature");
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
		addParam(createParamCentered<RoundHugeBlackKnob>(mm2px(Vec(23.761, 23.1)), module, Loom::COARSE_TUNE_KNOB_PARAM));

		// Regular knobs
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(46.643, 27.816)), module, Loom::HARM_COUNT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(61.994, 27.816)), module, Loom::HARM_DENSITY_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(77.345, 27.816)), module, Loom::HARM_HARM_STRIDE_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(92.695, 27.816)), module, Loom::HARM_SHIFT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(33.611, 61.274)), module, Loom::SPECTRAL_PIVOT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(53.975, 61.274)), module, Loom::SPECTRAL_TILT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(74.34, 61.274)), module, Loom::SPECTRAL_INTENSITY_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(33.611, 81.596)), module, Loom::HARMONIC_PIVOT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(53.975, 81.596)), module, Loom::HARMONIC_TILT_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(74.34, 81.627)), module, Loom::HARMONIC_INTENSITY_KNOB_PARAM));

		// Trimpots
		addParam(createParamCentered<Trimpot>(mm2px(Vec(46.643, 44.874)), module, Loom::HARM_COUNT_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(61.994, 44.874)), module, Loom::HARM_DENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(77.345, 44.874)), module, Loom::HARM_STRIDE_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(92.695, 44.874)), module, Loom::HARM_SHIFT_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(92.743, 61.724)), module, Loom::HARMONIC_INTENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(92.743, 81.627)), module, Loom::SPECTRAL_INTENSITY_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(7.817, 99.337)), module, Loom::PM_ATTENUVERTER_PARAM));
		addParam(createParamCentered<Trimpot>(mm2px(Vec(7.817, 42)), module, Loom::FM_ATTENUVERTER_PARAM));

		// Vertical switches
		addParam(createParamCentered<CKSS>(mm2px(Vec(6.385, 22.916)), module, Loom::RANGE_SWITCH_PARAM));
		addParam(createParamCentered<CKSSThree>(mm2px(Vec(13.822, 61.274)), module, Loom::HARMONIC_CURVE_SWITCH_PARAM));
		addParam(createParamCentered<CKSSThree>(mm2px(Vec(13.822, 81.177)), module, Loom::SPECTRAL_CURVE_SWITCH_PARAM));

		// Horizontal switches
		addParam(createParamCentered<CKSSHorizontal>(mm2px(Vec(48.192, 16.348)), module, Loom::CHARACTER_SWITCH_PARAM));
		addParam(createParamCentered<CKSSHorizontal>(mm2px(Vec(69.669, 16.348)), module, Loom::CONTINUOUS_STRIDE_SWITCH_PARAM));
		addParam(createParamCentered<CKSSHorizontal>(mm2px(Vec(91.147, 16.348)), module, Loom::INTERPOLATION_SWITCH_PARAM));
		addParam(createParamCentered<CKSSHorizontal>(mm2px(Vec(23.753, 42.057)), module, Loom::LIN_EXP_FM_SWITCH_PARAM));

		// Inputs
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(27.306, 99.133)), module, Loom::SYNC_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(37.051, 99.133)), module, Loom::HARM_COUNT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(46.795, 99.133)), module, Loom::HARM_STRIDE_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(56.54, 99.133)), module, Loom::SPECTRAL_PIVOT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(66.284, 99.133)), module, Loom::SPECTRAL_TILT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(76.029, 99.133)), module, Loom::SPECTRAL_INTENSITY_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(7.817, 111.779)), module, Loom::PM_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(17.562, 111.779)), module, Loom::FM_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(27.306, 111.779)), module, Loom::PITCH_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(37.051, 111.779)), module, Loom::HARM_DENSITY_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(46.795, 111.779)), module, Loom::HARM_SHIFT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(56.54, 111.779)), module, Loom::WAVESHAPE_PIVOT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(66.284, 111.779)), module, Loom::WAVESHAPE_TILT_CV_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(76.029, 111.779)), module, Loom::WAVESHAPE_INTENSITY_CV_INPUT));

		// Outputs
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(85.773, 99.133)), module, Loom::SQUARE_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(95.518, 99.133)), module, Loom::MAIN_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(85.773, 111.779)), module, Loom::SUB_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(95.518, 111.779)), module, Loom::QUADRATURE_OUT_OUTPUT));

		// Multi-colored LEDs
		addChild(createLightCentered<MediumLight<GreenRedLight>>(mm2px(Vec(33.4, 13.3)), module, Loom::OSCILLATOR_LED_LIGHT));

		// Single color LEDs
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(69.75, 25.882)), module, Loom::STRIDE_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(82.201, 53.808)), module, Loom::S_LED_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(86.259, 53.808)), module, Loom::S_LED_2_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(90.317, 53.808)), module, Loom::S_LED_3_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(94.376, 53.808)), module, Loom::S_LED_4_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(98.434, 53.808)), module, Loom::S_LED_5_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(82.064, 73.834)), module, Loom::H_LED_1_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(86.122, 73.834)), module, Loom::H_LED_2_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(90.181, 73.834)), module, Loom::H_LED_3_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(94.239, 73.834)), module, Loom::H_LED_4_LIGHT));
		addChild(createLightCentered<SmallSimpleLight<BlueLight>>(mm2px(Vec(98.297, 73.834)), module, Loom::H_LED_5_LIGHT));
	}
};


Model* modelLoom = createModel<Loom, LoomWidget>("RigatoniModular-Loom");