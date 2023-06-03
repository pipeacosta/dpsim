/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/DP/DP_Ph3_SynchronGeneratorDQ.h>

namespace CPS {
namespace DP {
namespace Ph3 {
	class SynchronGeneratorDQTrapez :
		public SynchronGeneratorDQ,
		public SharedFactory<SynchronGeneratorDQTrapez> {
	public:
		SynchronGeneratorDQTrapez(String uid, String name, Logger::Level loglevel = Logger::Level::off);
		SynchronGeneratorDQTrapez(String name, Logger::Level loglevel = Logger::Level::off);

		void setMultisamplingRate(Int rate);

		// #### MNA Section ####
		void mnaCompInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector);

		/// Add MNA pre step dependencies
		void mnaCompAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) override;
		void mnaCompPreStep(Real time, Int timeStepCount) override;

	protected:
		///
		Int mMultisamplingRate = 1;

		// #### Trapezoidal Section ####

		/// Performs an Euler forward step with the state space model of a synchronous generator
		/// to calculate the flux and current from the voltage vector in per unit.
		void stepInPerUnit(Real time);
	};
}
}
}
