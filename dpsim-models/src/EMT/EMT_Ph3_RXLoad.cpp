/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_RXLoad.h>

using namespace CPS;

EMT::Ph3::RXLoad::RXLoad(String uid, String name, Logger::Level logLevel)
	: CompositePowerComp<Real>(uid, name, true, true, logLevel),
	mActivePower(Attribute<Matrix>::create("P", mAttributes)),
	mReactivePower(Attribute<Matrix>::create("Q", mAttributes)),
	mNomVoltage(Attribute<Real>::create("V_nom", mAttributes)) {
	mPhaseType = PhaseType::ABC;
	setTerminalNumber(1);

	mSLog->info("Create {} {}", this->type(), name);
	**mIntfVoltage = Matrix::Zero(3, 1);
	**mIntfCurrent = Matrix::Zero(3, 1);
	mSLog->flush();
}

EMT::Ph3::RXLoad::RXLoad(String name,
	Logger::Level logLevel)
	: RXLoad(name, name, logLevel) {
}

EMT::Ph3::RXLoad::RXLoad(String name,
	Matrix activePower, Matrix reactivePower, Real volt,
	Logger::Level logLevel)
	: RXLoad(name, logLevel) {
	**mActivePower = activePower;
	**mReactivePower = reactivePower;
	mPower = MatrixComp::Zero(3,3);
	mPower <<
		Complex((**mActivePower)(0, 0), (**mReactivePower)(0, 0)), Complex((**mActivePower)(0, 1), (**mReactivePower)(0, 1)), Complex((**mActivePower)(0, 2), (**mReactivePower)(0, 2)),
		Complex((**mActivePower)(1, 0), (**mReactivePower)(1, 0)), Complex((**mActivePower)(1, 1), (**mReactivePower)(1, 1)), Complex((**mActivePower)(1, 2), (**mReactivePower)(1, 2)),
		Complex((**mActivePower)(2, 0), (**mReactivePower)(2, 0)), Complex((**mActivePower)(2, 1), (**mReactivePower)(2, 1)), Complex((**mActivePower)(2, 2), (**mReactivePower)(2, 2));

	**mNomVoltage = volt;
	initPowerFromTerminal = false;
}

void EMT::Ph3::RXLoad::setParameters(Matrix activePower, Matrix reactivePower, Real volt) {
	**mActivePower = activePower;
	**mReactivePower = reactivePower;

	// complex power
	mPower = MatrixComp::Zero(3, 3);
	mPower(0, 0) = { (**mActivePower)(0, 0), (**mReactivePower)(0, 0) };
	mPower(1, 1) = { (**mActivePower)(1, 1), (**mReactivePower)(1, 1) };
	mPower(2, 2) = { (**mActivePower)(2, 2), (**mReactivePower)(2, 2) };

	**mNomVoltage = volt;

	mSLog->info("\nActive Power [W]: {}"
			"\nReactive Power [VAr]: {}",
			Logger::matrixToString(**mActivePower),
			Logger::matrixToString(**mReactivePower));
	mSLog->info("Nominal Voltage={} [V]", **mNomVoltage);

	initPowerFromTerminal = false;
}

SimPowerComp<Real>::Ptr EMT::Ph3::RXLoad::clone(String name) {
	// everything set by initializeFromNodesAndTerminals
	return RXLoad::make(name, mLogLevel);
}

void EMT::Ph3::RXLoad::initializeFromNodesAndTerminals(Real frequency) {

		if (initPowerFromTerminal) {
		**mActivePower = Matrix::Zero(3, 3);
		(**mActivePower)(0, 0) = mTerminals[0]->singleActivePower() / 3.;
		(**mActivePower)(1, 1) = mTerminals[0]->singleActivePower() / 3.;
		(**mActivePower)(2, 2) = mTerminals[0]->singleActivePower() / 3.;

		**mReactivePower = Matrix::Zero(3, 3);
		(**mReactivePower)(0, 0) = mTerminals[0]->singleReactivePower() / 3.;
		(**mReactivePower)(1, 1) = mTerminals[0]->singleReactivePower() / 3.;
		(**mReactivePower)(2, 2) = mTerminals[0]->singleReactivePower() / 3.;

		// complex power
		mPower = MatrixComp::Zero(3, 3);
		mPower(0, 0) = { (**mActivePower)(0, 0), (**mReactivePower)(0, 0) };
		mPower(1, 1) = { (**mActivePower)(1, 1), (**mReactivePower)(1, 1) };
		mPower(2, 2) = { (**mActivePower)(2, 2), (**mReactivePower)(2, 2) };

		**mNomVoltage = std::abs(mTerminals[0]->initialSingleVoltage());

		mSLog->info("\nActive Power [W]: {}"
					"\nReactive Power [VAr]: {}",
					Logger::matrixToString(**mActivePower),
					Logger::matrixToString(**mReactivePower));
		mSLog->info("Nominal Voltage={} [V]", **mNomVoltage);
	}

	if ((**mActivePower)(0,0) != 0) {
		mResistance = std::pow(**mNomVoltage / sqrt(3), 2) * (**mActivePower).inverse();
		mSubResistor = std::make_shared<EMT::Ph3::Resistor>(**mName + "_res", mLogLevel);
		mSubResistor->setParameters(mResistance);
		mSubResistor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubResistor->initialize(mFrequencies);
		mSubResistor->initializeFromNodesAndTerminals(frequency);
		addMNASubComponent(mSubResistor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
	}

	if ((**mReactivePower)(0, 0) != 0)
		mReactance = std::pow(**mNomVoltage / sqrt(3), 2) * (**mReactivePower).inverse();
	else
		mReactance = Matrix::Zero(1, 1);

	if (mReactance(0,0) > 0) {
		mInductance = mReactance / (2 * PI * frequency);

		mSubInductor = std::make_shared<EMT::Ph3::Inductor>(**mName + "_ind", mLogLevel);
		mSubInductor->setParameters(mInductance);
		mSubInductor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubInductor->initialize(mFrequencies);
		mSubInductor->initializeFromNodesAndTerminals(frequency);
		addMNASubComponent(mSubInductor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
	}
	else if (mReactance(0,0) < 0) {
		mCapacitance = -1 / (2 * PI * frequency) * mReactance.inverse();

		mSubCapacitor = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_cap", mLogLevel);
		mSubCapacitor->setParameters(mCapacitance);
		mSubCapacitor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubCapacitor->initialize(mFrequencies);
		mSubCapacitor->initializeFromNodesAndTerminals(frequency);
		addMNASubComponent(mSubCapacitor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
	}

	MatrixComp vInitABC = MatrixComp::Zero(3, 1);
	vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * mTerminals[0]->initialSingleVoltage();
	vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
	vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
	**mIntfVoltage = vInitABC.real();

	MatrixComp iInitABC = MatrixComp::Zero(3, 1);
	// v i^T* = S
	// v^T v i^T* = v^T S
	// i^T*= (|v|^2)^(-1) v^T S

	Complex v_ = vInitABC(0, 0)*vInitABC(0, 0) + vInitABC(1, 0)*vInitABC(1, 0) + vInitABC(2, 0)*vInitABC(2, 0);
	MatrixComp rhs_ = Complex(1, 0) / v_ * vInitABC.transpose() * mPower;
	iInitABC = rhs_.conjugate().transpose();
	**mIntfCurrent = iInitABC.real();

	mSLog->info(
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nActive Power: {:s}"
		"\nReactive Power: {:s}"
		"\nResistance: {:s}"
		"\nReactance: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::matrixToString(**mIntfVoltage),
		Logger::matrixToString(**mIntfCurrent),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
		Logger::matrixToString(**mActivePower),
		Logger::matrixToString(**mReactivePower),
		Logger::matrixToString(mResistance),
		Logger::matrixToString(mReactance));
	mSLog->flush();

}

void EMT::Ph3::RXLoad::mnaParentAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
	modifiedAttributes.push_back(mRightVector);
};

void EMT::Ph3::RXLoad::mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(mIntfCurrent);
	modifiedAttributes.push_back(mIntfVoltage);
};

void EMT::Ph3::RXLoad::mnaParentPreStep(Real time, Int timeStepCount) {
	mnaApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::RXLoad::mnaParentPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	mnaUpdateVoltage(**leftVector);
	mnaUpdateCurrent(**leftVector);
}

void EMT::Ph3::RXLoad::mnaUpdateVoltage(const Matrix& leftVector) {
	**mIntfVoltage = Matrix::Zero(3, 1);
	(**mIntfVoltage)(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
	(**mIntfVoltage)(1, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
	(**mIntfVoltage)(2, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
}

void EMT::Ph3::RXLoad::mnaUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = Matrix::Zero(3, 1);
	for (auto& subc : mSubComponents) {
		**mIntfCurrent += subc->intfCurrent();
	}
}
