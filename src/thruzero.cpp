#include "plugin.hpp"


struct ThruZero : Module {
	enum ParamId {
		THRESHOLD_KNOB_PARAM,
		IN_LVL_KNOB_PARAM,
		PARAMS_LEN
	};
	enum InputId {
		IN_LVL_IN_INPUT,
		THRESH_IN_INPUT,
		CV_IN_INPUT,
		RAMP_IN_INPUT,
		INPUTS_LEN
	};
	enum OutputId {
		FWD_OUT_OUTPUT,
		CV_OUT_OUTPUT,
		SINE_OUT_OUTPUT,
		RAMP_OUT_OUTPUT,
		OUTPUTS_LEN
	};
	enum LightId {
		FWD_LED_LIGHT,
		BACK_LED_LIGHT,
		LIGHTS_LEN
	};

	ThruZero() {
		config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);
		configParam(THRESHOLD_KNOB_PARAM, -5.f, 5.f, 0.f, "Thru-zero threshold");
		configParam(IN_LVL_KNOB_PARAM, -2.f, 2.f, 1.f, "Input ramp level");
		configInput(IN_LVL_IN_INPUT, "Input level CV");
		configInput(THRESH_IN_INPUT, "Threshold CV");
		configInput(CV_IN_INPUT, "1V/oct CV input");
		configInput(RAMP_IN_INPUT, "Ramp input");
		configOutput(FWD_OUT_OUTPUT, "Forward direction gate");
		configOutput(CV_OUT_OUTPUT, "1V/oct CV output");
		configOutput(SINE_OUT_OUTPUT, "Sine output");
		configOutput(RAMP_OUT_OUTPUT, "Ramp output");

		lightDivider.setDivision(32);
	}

	dsp::RCFilter spikeLPF;
	dsp::RCFilter spikeIsolatorHPF;
	float spikeEnvF = 0.f;
	float sampledCompThresh = 5.f;
	float lastRampVal = 0.f;
	bool lastFwd = true;
	dsp::ClockDivider lightDivider;

	void process(const ProcessArgs& args) override {
		float thresh = clamp(params[THRESHOLD_KNOB_PARAM].getValue() + inputs[THRESH_IN_INPUT].getVoltage(), -5.f, 5.0f);
		float inLvl = clamp(params[IN_LVL_KNOB_PARAM].getValue() + .4f * inputs[IN_LVL_IN_INPUT].getVoltage(), -2.f, 2.f);

		// Calculate CV out and fwd/back
		float offsetCV = inputs[CV_IN_INPUT].getVoltage() - thresh;
		bool newFwd = offsetCV >= 0.f;
		offsetCV = std::fabs(offsetCV);
		float ledBrightness = std::min(offsetCV * .2f, 1.f);
		float outCV = offsetCV + thresh;

		// Calculate PM ramp output
		float scaledRamp = inputs[RAMP_IN_INPUT].getVoltage() * inLvl * (lastFwd ? 1.f : -1.f);
		// if (scaledRamp > .5f) scaledRamp -= 10.f;
		// if (scaledRamp < -.5f) scaledRamp += 10.f;
		float negComp = (scaledRamp > sampledCompThresh) ? -10.f : 0.f;
		float posComp = negComp + 10.f;
		float rampOut = clamp(scaledRamp + crossfade(posComp, negComp, .1f * (sampledCompThresh + 5.f)), -5.f, 5.f);

		// Remove spike from phase-shifted ramp
		spikeIsolatorHPF.setCutoff(10000 / args.sampleRate);
		spikeIsolatorHPF.process(scaledRamp);
		float hpfRect = std::fabs(spikeIsolatorHPF.highpass());
		spikeEnvF = (hpfRect >= spikeEnvF) ? (hpfRect * .8f + spikeEnvF * .2f) : spikeEnvF * .5f;
		spikeLPF.setCutoff(((spikeEnvF >= 2.5f) ? 500 : 20000) / args.sampleRate);
		spikeLPF.process(rampOut);
		// rampOut = spikeLPF.lowpass();

		// Check for FM polarity crossing and change crossing sample
		if (lastFwd != newFwd) {
			float crossingSample = rampOut;
			float unshiftedInvertedSample = -scaledRamp;
			float shiftAmt = unshiftedInvertedSample - crossingSample;
			sampledCompThresh = shiftAmt - 5.f;
			if (sampledCompThresh > 5.f) sampledCompThresh -= 10.f;
			if (sampledCompThresh < -5.f) sampledCompThresh += 10.f;
		}

		// Outputs
		outputs[FWD_OUT_OUTPUT].setVoltage(newFwd ? 10.f : 0.f);
		outputs[CV_OUT_OUTPUT].setVoltage(outCV);
		outputs[RAMP_OUT_OUTPUT].setVoltage(rampOut);
		outputs[SINE_OUT_OUTPUT].setVoltage(5.f * std::sin(2 * float(M_PI) * (.1f * (rampOut + 5.f))));

		// Update internal state
		lastFwd = newFwd;
		lastRampVal = scaledRamp;

		if (lightDivider.process()) {
			float lightTime = args.sampleTime * lightDivider.getDivision();
			lights[FWD_LED_LIGHT].setBrightnessSmooth(newFwd ? ledBrightness : 0.f, lightTime);
			lights[BACK_LED_LIGHT].setBrightnessSmooth(newFwd ? 0.f : ledBrightness, lightTime);
		}
	}
};


struct ThruZeroWidget : ModuleWidget {
	ThruZeroWidget(ThruZero* module) {
		setModule(module);
		setPanel(createPanel(asset::plugin(pluginInstance, "res/ThruZero.svg")));

		addChild(createWidget<ScrewSilver>(Vec(0, 0)));
		addChild(createWidget<ScrewSilver>(Vec(0, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		addParam(createParamCentered<RoundBigBlackKnob>(mm2px(Vec(9.985, 32.445)), module, ThruZero::THRESHOLD_KNOB_PARAM));
		addParam(createParamCentered<RoundBigBlackKnob>(mm2px(Vec(10.16, 60.304)), module, ThruZero::IN_LVL_KNOB_PARAM));

		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.058, 81.75)), module, ThruZero::IN_LVL_IN_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.058, 91.75)), module, ThruZero::THRESH_IN_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.058, 101.75)), module, ThruZero::CV_IN_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.058, 111.75)), module, ThruZero::RAMP_IN_INPUT));

		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(14.802, 81.75)), module, ThruZero::FWD_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(14.802, 91.75)), module, ThruZero::CV_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(14.802, 101.75)), module, ThruZero::SINE_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(14.802, 111.75)), module, ThruZero::RAMP_OUT_OUTPUT));

		addChild(createLightCentered<MediumLight<GreenLight>>(mm2px(Vec(8.957, 12.194)), module, ThruZero::FWD_LED_LIGHT));
		addChild(createLightCentered<MediumLight<RedLight>>(mm2px(Vec(8.957, 16.139)), module, ThruZero::BACK_LED_LIGHT));
	}
};


Model* modelThruZero = createModel<ThruZero, ThruZeroWidget>("ThruZero");