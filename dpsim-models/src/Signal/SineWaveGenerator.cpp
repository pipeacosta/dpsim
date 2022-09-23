/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/Signal/SineWaveGenerator.h>

using namespace CPS;

Signal::SineWaveGenerator::SineWaveGenerator(String name, Logger::Level logLevel)
			: SignalGenerator(name, logLevel),
			mMagnitude(Attribute<Real>::createDynamic("V_max", mAttributes)) { }

void Signal::SineWaveGenerator::setParameters(Complex initialPhasor, Real frequency /*= -1*/) {
    **mMagnitude = Math::abs(initialPhasor);
    mInitialPhase = Math::phase(initialPhasor);
	**mFreq = frequency;
	**mSigOut = initialPhasor;

	mSLog->info(	"Parameters: \n"
					 "Sine wave magnitude: {} [V] \n"
					 "Sine wave initial phase: {} [rad] \n"
					 "Sine wave frequency: {} [Hz] \n",
					**mMagnitude, mInitialPhase, **mFreq);				 
}

void Signal::SineWaveGenerator::step(Real time) {
	**mSigOut = Complex(
		**mMagnitude * cos(time * 2.*PI* **mFreq + mInitialPhase),
		**mMagnitude * sin(time * 2.*PI* **mFreq + mInitialPhase));
}