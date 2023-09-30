#include "plugin.hpp"

struct SpikeProcessor {
	dsp::RCFilter spikeLPF;
	dsp::RCFilter spikeIsolatorHPF;
	float spikeEnvF = 0.f;

	SpikeProcessor() {}

	float process(float sampleRate, float ramp, float derived) {
		spikeIsolatorHPF.setCutoff(10000 / sampleRate);
		spikeIsolatorHPF.process(ramp);
		float hpfRect = std::fabs(spikeIsolatorHPF.highpass());
		spikeEnvF = (hpfRect >= spikeEnvF) ? (hpfRect * .75f + spikeEnvF * .25f) : spikeEnvF * .25f;
		spikeLPF.setCutoff(((spikeEnvF >= 2.5f) ? 500 : 20000) / sampleRate);
		spikeLPF.process(derived);
		return spikeLPF.lowpass();
	}
};

struct ThruZero : Module {
	enum ParamId {
		THRESHOLD_KNOB_PARAM,
		IN_LVL_KNOB_PARAM,
		SMOOTH_SWITCH_PARAM,
		FLAVOR_SWITCH_PARAM,
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
		configSwitch(SMOOTH_SWITCH_PARAM, 0.f, 1.f, 1.f, "Enable spike smoothing", {"Off", "On"});
		configSwitch(FLAVOR_SWITCH_PARAM, 0.f, 1.f, 0.f, "Enable classic XWY flavor", {"Off", "On"});
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

	SpikeProcessor rampSpikeProcessor;
	SpikeProcessor sinSpikeProcessor;
	float sampledCompThresh = 5.f;
	bool lastFwd = true;
	dsp::ClockDivider lightDivider;

	void process(const ProcessArgs& args) override {
		float thresh = clamp(params[THRESHOLD_KNOB_PARAM].getValue() + inputs[THRESH_IN_INPUT].getVoltage(), -5.f, 5.0f);
		float inLvl = clamp(params[IN_LVL_KNOB_PARAM].getValue() + .4f * inputs[IN_LVL_IN_INPUT].getVoltage(), -2.f, 2.f);
		bool xwyFlavor = params[FLAVOR_SWITCH_PARAM].getValue() > .5f;
		bool smoothing = params[SMOOTH_SWITCH_PARAM].getValue() > .5f;

		// Calculate modulator CV output from modulator CV input based on "thru-zero" threshold
		float offsetCV = inputs[CV_IN_INPUT].getVoltage() - thresh;
		bool newFwd = offsetCV >= 0.f;
		offsetCV = std::fabs(offsetCV);
		float outCV = offsetCV + thresh;

		// LED display on panel will show whether we're outputting a positive or negative modulation signal
		float ledBrightness = std::min(offsetCV * .2f, 1.f);

		// Calculate and remove spike from phase-shifted ramp output, which could be
		// flipped if we've done an odd # of "thru-zero" crossings
		float scaledRamp = inputs[RAMP_IN_INPUT].getVoltage() * inLvl * (lastFwd ? 1.f : -1.f);
		float rampOut = phaseShiftRamp(scaledRamp, sampledCompThresh);
		if (smoothing) rampOut = rampSpikeProcessor.process(args.sampleRate, scaledRamp, rampOut);

		// Check for FM polarity crossing and change crossing sample
		if (lastFwd != newFwd) {
			// Do some demonic math that essentially finds the right phase shift for the next
			// inversion of our input ramp. We're sample-and-hold-ing a comparator threshold
			// here rather than saving a phase shift since the comparator threshold is more
			// immediate to work with
			float crossingSample = rampOut;
			float unshiftedInvertedSample = -scaledRamp;
			float shiftAmt = unshiftedInvertedSample - crossingSample;
			sampledCompThresh = shiftAmt - 5.f;
			// As much as I tried to simplify the algebra above, annoyingly the range still goes
			// outside of -5V to 5V due to the analytical methods I used, so here's a nice hack
			if (sampledCompThresh > 5.f) sampledCompThresh -= 10.f;
			if (sampledCompThresh < -5.f) sampledCompThresh += 10.f;
		}

		// Do the X Without Y version
		if (xwyFlavor) rampOut = newFwd ? rampOut : -rampOut;

		// Straightforward outputs
		outputs[FWD_OUT_OUTPUT].setVoltage(newFwd ? 10.f : 0.f);
		outputs[CV_OUT_OUTPUT].setVoltage(outCV);
		outputs[RAMP_OUT_OUTPUT].setVoltage(rampOut);

		// Was too lazy to implement a better saw-to-sin shaper so here's a quadratic approximation
		float sin = rampOut * 3.f;
		float absSin = std::fabs(sin);
		if (absSin > 5.f) sin -= .15f * std::pow(absSin - 5.f, 2.f) * (sin > 0.f ? 1.f : -1.f);
		if (smoothing) sin = sinSpikeProcessor.process(args.sampleRate, rampOut, sin);
		outputs[SINE_OUT_OUTPUT].setVoltage(.75f * sin);

		// Update internal state
		lastFwd = newFwd;

		if (lightDivider.process()) {
			float lightTime = args.sampleTime * lightDivider.getDivision();
			lights[FWD_LED_LIGHT].setBrightnessSmooth(newFwd ? ledBrightness : 0.f, lightTime);
			lights[BACK_LED_LIGHT].setBrightnessSmooth(newFwd ? 0.f : ledBrightness, lightTime);
		}
	}

	// My favorite way to patch up a ramp phase shifter on my modular rig,
	// using some gain, offsets, a comparator, and a crossfader
	float phaseShiftRamp(float sample, float thresh) {
		float negComp = (sample > thresh) ? -10.f : 0.f;
		float posComp = negComp + 10.f;
		return clamp(sample + crossfade(posComp, negComp, .1f * (thresh + 5.f)), -5.f, 5.f);
	}
};

struct ThruZeroWidget : ModuleWidget {
	ThruZeroWidget(ThruZero* module) {
		setModule(module);
		setPanel(createPanel(asset::plugin(pluginInstance, "res/panels/thruzero.svg")));

		addChild(createWidget<ScrewSilver>(Vec(0, 0)));
		addChild(createWidget<ScrewSilver>(Vec(0, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		addParam(createParamCentered<RoundBigBlackKnob>(mm2px(Vec(10.16, 40.911)), module, ThruZero::THRESHOLD_KNOB_PARAM));
		addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(10.16, 64.403)), module, ThruZero::IN_LVL_KNOB_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(6.22, 16.121)), module, ThruZero::SMOOTH_SWITCH_PARAM));
		addParam(createParamCentered<CKSS>(mm2px(Vec(14.279, 16.121)), module, ThruZero::FLAVOR_SWITCH_PARAM));

		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.058, 81.75)), module, ThruZero::IN_LVL_IN_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.058, 91.75)), module, ThruZero::THRESH_IN_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.058, 101.75)), module, ThruZero::CV_IN_INPUT));
		addInput(createInputCentered<PJ301MPort>(mm2px(Vec(5.058, 111.75)), module, ThruZero::RAMP_IN_INPUT));

		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(14.802, 81.75)), module, ThruZero::FWD_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(14.802, 91.75)), module, ThruZero::CV_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(14.802, 101.75)), module, ThruZero::SINE_OUT_OUTPUT));
		addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(14.802, 111.75)), module, ThruZero::RAMP_OUT_OUTPUT));

		addChild(createLightCentered<MediumLight<GreenLight>>(mm2px(Vec(14.129, 27.780)), module, ThruZero::FWD_LED_LIGHT));
		addChild(createLightCentered<MediumLight<RedLight>>(mm2px(Vec(6.191, 27.780)), module, ThruZero::BACK_LED_LIGHT));
	}
};


Model* modelThruZero = createModel<ThruZero, ThruZeroWidget>("RigatoniModular-ThruZero");