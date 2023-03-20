/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/Solver/PFSolverInterfaceBranch.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/SP/SP_Ph1_Resistor.h>
#include <dpsim-models/SP/SP_Ph1_Inductor.h>
#include <dpsim-models/SP/SP_Ph1_Capacitor.h>
#include <dpsim-models/Base/Base_Ph1_Transformer.h>

namespace CPS {
namespace SP {
namespace Ph1 {
	/// Transformer that includes an inductance and resistance
	class Transformer :
		public CompositePowerComp<Complex>,
		public Base::Ph1::Transformer,
		public SharedFactory<Transformer>,
		public PFSolverInterfaceBranch {

	private:
		/// Internal resistor to model losses
		std::shared_ptr<SP::Ph1::Resistor> mSubResistor;
		/// Internal inductor to model losses
		std::shared_ptr<SP::Ph1::Inductor> mSubInductor;

		/// Internal parallel resistance 1 as snubber
		std::shared_ptr<SP::Ph1::Resistor> mSubSnubResistor1;
		/// Internal parallel resistance 2 as snubber
		std::shared_ptr<SP::Ph1::Resistor> mSubSnubResistor2;
		/// Internal parallel capacitance 1 as snubber
		std::shared_ptr<SP::Ph1::Capacitor> mSubSnubCapacitor1;
		/// Internal parallel capacitance 2 as snubber
		std::shared_ptr<SP::Ph1::Capacitor> mSubSnubCapacitor2;

		/// Snubber resistance 1 [Ohm]
		Real mSnubberResistance1;
		/// Snubber resistance 2 [Ohm]
		Real mSnubberResistance2;
		/// Snubber capacitance 1 [F]
		Real mSnubberCapacitance1;
		/// Snubber capacitance 2 [F]
		Real mSnubberCapacitance2;

        /// Transformer ratio magnitude
		Real mRatioAbs = 1;
        /// Transformer ratio pase [deg]
		Real mRatioPhase = 0;
		/// Nominal omega
		Real mNominalOmega;
		/// Reactance [Ohm]
		Real mReactance;

		/// Leakage
		Complex mLeakage;

		/// base apparent power[VA]
		Real mBaseApparentPower;
		/// base impedance [ohm]
		Real mBaseImpedance;
        /// base inductance [H]
        Real mBaseInductance;
		/// base admittance [S]
		Real mBaseAdmittance;
		///base omega [1/s]
		Real mBaseOmega;

		///base current [A]
		Real mBaseCurrent;

		/// resistance
		Real mResistancePerUnit;
		/// reactance
		Real mReactancePerUnit;
        /// inductance
        Real mInductancePerUnit;
		/// leakage impedance
		Complex mLeakagePerUnit;
		/// magnetizing impedance
		Complex mMagnetizingPerUnit;
        /// transformer ratio
        Real mRatioAbsPerUnit;

		// #### Admittance matrix stamp ####
		MatrixComp mY_element;

		/// Boolean for considering resistive losses with sub resistor
		Bool mWithResistiveLosses;
	public:
		/// base voltage [V]
		const Attribute<Real>::Ptr mBaseVoltage;

		// #### Power flow results ####
		/// branch Current flow [A], coef(0) has data from node 0, coef(1) from node 1.
		const Attribute<MatrixComp>::Ptr mCurrent;

		/// branch active powerflow [W], coef(0) has data from node 0, coef(1) from node 1.
		const Attribute<Matrix>::Ptr mActivePowerBranch;

		/// branch reactive powerflow [Var], coef(0) has data from node 0, coef(1) from node 1.
		const Attribute<Matrix>::Ptr mReactivePowerBranch;
		/// nodal active power injection
		const Attribute<Real>::Ptr mActivePowerInjection;
		/// nodal reactive power injection
		const Attribute<Real>::Ptr mReactivePowerInjection;

		/// Defines UID, name and logging level
		Transformer(String uid, String name,
			Logger::Level logLevel = Logger::Level::off, Bool withResistiveLosses = false);
		/// Defines name and logging level
		Transformer(String name, Logger::Level logLevel = Logger::Level::off)
			: Transformer(name, name, logLevel) { }

		SimPowerComp<Complex>::Ptr clone(String name) override;

		// #### General ####
		/// Set transformer specific parameters (without rated power)
		void setParameters(Real nomVoltageEnd1, Real nomVoltageEnd2, Real ratioAbs, Real ratioPhase, Real resistance, Real inductance);
		/// Set transformer specific parameters
		void setParameters(Real nomVoltageEnd1, Real nomVoltageEnd2, Real ratedPower, Real ratioAbs, Real ratioPhase, Real resistance, Real inductance);
		/// Initializes component from power flow data
		void initializeFromNodesAndTerminals(Real frequency) override;

		// #### Powerflow section ####
		/// Set base voltage
		void setBaseVoltage(Real baseVoltage);
		/// Initializes component from power flow data
		void calculatePerUnitParameters(Real baseApparentPower, Real baseOmega);
		/// Stamps admittance matrix
		void pfApplyAdmittanceMatrixStamp(SparseMatrixCompRow & Y) override;
		/// updates branch current and power flow, input pu value, update with real value
		void updateBranchFlow(VectorComp& current, VectorComp& powerflow);
		/// stores nodal injection power in this line object
		void storeNodalInjection(Complex powerInjection);

		// #### Getter ####
		/// get admittance matrix
		MatrixComp Y_element();

		// #### MNA Section ####
		/// Initializes internal variables of the component
		void mnaParentInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) override;
		/// Stamps system matrix
		void mnaCompApplySystemMatrixStamp(SparseMatrixRow& systemMatrix) override;
		/// Updates internal current variable of the component
		void mnaCompUpdateCurrent(const Matrix& leftVector) override;
		/// Updates internal voltage variable of the component
		void mnaCompUpdateVoltage(const Matrix& leftVector) override;
		/// MNA pre step operations
		void mnaParentPreStep(Real time, Int timeStepCount) override;
		/// MNA post step operations
		void mnaParentPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) override;
		/// Add MNA pre step dependencies
		void mnaParentAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) override;
		/// Add MNA post step dependencies
		void mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) override;

    };
}
}
}
