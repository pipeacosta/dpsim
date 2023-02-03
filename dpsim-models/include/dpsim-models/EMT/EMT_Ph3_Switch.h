/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/SimPowerComp.h>
#include <dpsim-models/Solver/MNASwitchInterface.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/Logger.h>
#include <dpsim-models/Base/Base_Ph3_Switch.h>

namespace CPS {
namespace EMT {
namespace Ph3 {
	/// \brief Dynamic phasor switch
	///
	/// The switch can be opened and closed.
	/// Each state has a specific resistance value.
	class Switch :
		public Base::Ph3::Switch,
		public SimPowerComp<Real>,
		public SharedFactory<Switch>,
		public MNASwitchInterface,
		public MNAInterface {
	protected:
	public:
		/// Defines UID, name, component parameters and logging level
		Switch(String uid, String name,	Logger::Level loglevel = Logger::Level::off);
		/// Defines name, component parameters and logging level
		Switch(String name, Logger::Level logLevel = Logger::Level::off)
			: Switch(name, name, logLevel) { }

		SimPowerComp<Real>::Ptr clone(String name);

		// #### General ####
		/// Initializes component from power flow data
		void initializeFromNodesAndTerminals(Real frequency);

		// #### General MNA section ####
		void mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector);
		/// Stamps system matrix
		void mnaApplySystemMatrixStamp(Matrix& systemMatrix);
		/// Stamps right side (source) vector
		void mnaApplyRightSideVectorStamp(Matrix& rightVector);
		/// Update interface voltage from MNA system result
		void mnaUpdateVoltage(const Matrix& leftVector);
		/// Update interface current from MNA system result
		void mnaUpdateCurrent(const Matrix& leftVector);

		// #### MNA section for switches ####
		/// Check if switch is closed
		Bool mnaIsClosed() override;
		/// Stamps system matrix considering the defined switch position
		void mnaApplySwitchSystemMatrixStamp(Bool closed, Matrix& systemMatrix, Int freqIdx);

		class MnaPostStep : public Task {
		public:
			MnaPostStep(Switch& switchRef, Attribute<Matrix>::Ptr leftSideVector) :
				Task(**switchRef.mName + ".MnaPostStep"), mSwitch(switchRef), mLeftVector(leftSideVector) {
				mAttributeDependencies.push_back(mLeftVector);
				mModifiedAttributes.push_back(mSwitch.mIntfVoltage);
				mModifiedAttributes.push_back(mSwitch.mIntfCurrent);
			}

			void execute(Real time, Int timeStepCount);

		private:
			Switch& mSwitch;
			Attribute<Matrix>::Ptr mLeftVector;
		};
	};
}
}
}
