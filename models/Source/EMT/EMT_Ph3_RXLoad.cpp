/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/EMT/EMT_Ph3_RXLoad.h>

using namespace CPS;

EMT::Ph3::RXLoad::RXLoad(String uid, String name,
	Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	mPhaseType = PhaseType::ABC;
	setTerminalNumber(1);

	mSLog->info("Create {} {}", this->type(), name);
	mIntfVoltage = Matrix::Zero(3, 1);
	mIntfCurrent = Matrix::Zero(3, 1);

	addAttribute<Matrix>("P", &mActivePower, Flags::read | Flags::write);
	addAttribute<Matrix>("Q", &mReactivePower, Flags::read | Flags::write);
	addAttribute<Real>("V_nom", &mNomVoltage, Flags::read | Flags::write);
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
	mActivePower = activePower;
	mReactivePower = reactivePower;
	mPower = MatrixComp::Zero(3,3);
	mPower <<
		Complex(mActivePower(0, 0), mReactivePower(0, 0)), Complex(mActivePower(0, 1), mReactivePower(0, 1)), Complex(mActivePower(0, 2), mReactivePower(0, 2)),
		Complex(mActivePower(1, 0), mReactivePower(1, 0)), Complex(mActivePower(1, 1), mReactivePower(1, 1)), Complex(mActivePower(1, 2), mReactivePower(1, 2)),
		Complex(mActivePower(2, 0), mReactivePower(2, 0)), Complex(mActivePower(2, 1), mReactivePower(2, 1)), Complex(mActivePower(2, 2), mReactivePower(2, 2));

	mNomVoltage = volt;
	initPowerFromTerminal = false;
}

void EMT::Ph3::RXLoad::setParameters(Matrix activePower, Matrix reactivePower, Real volt) {
	mActivePower = activePower;
	mReactivePower = reactivePower;

	// complex power
	mPower = MatrixComp::Zero(3, 3);
	mPower(0, 0) = { mActivePower(0, 0), mReactivePower(0, 0) };
	mPower(1, 1) = { mActivePower(1, 1), mReactivePower(1, 1) };
	mPower(2, 2) = { mActivePower(2, 2), mReactivePower(2, 2) };

	mNomVoltage = volt;

	mSLog->info("\nActive Power [W]: {}"
			"\nReactive Power [VAr]: {}",
			Logger::matrixToString(mActivePower),
			Logger::matrixToString(mReactivePower));
	mSLog->info("Nominal Voltage={} [V]", mNomVoltage);

	initPowerFromTerminal = false;
}

SimPowerComp<Real>::Ptr EMT::Ph3::RXLoad::clone(String name) {
	// everything set by initializeFromNodesAndTerminals
	return RXLoad::make(name, mLogLevel);
}

void EMT::Ph3::RXLoad::initializeFromNodesAndTerminals(Real frequency) {

		if (initPowerFromTerminal) {
		mActivePower = Matrix::Zero(3, 3);
		mActivePower(0, 0) = mTerminals[0]->singleActivePower() / 3.;
		mActivePower(1, 1) = mTerminals[0]->singleActivePower() / 3.;
		mActivePower(2, 2) = mTerminals[0]->singleActivePower() / 3.;

		mReactivePower = Matrix::Zero(3, 3);
		mReactivePower(0, 0) = mTerminals[0]->singleReactivePower() / 3.;
		mReactivePower(1, 1) = mTerminals[0]->singleReactivePower() / 3.;
		mReactivePower(2, 2) = mTerminals[0]->singleReactivePower() / 3.;

		// complex power
		mPower = MatrixComp::Zero(3, 3);
		mPower(0, 0) = { mActivePower(0, 0), mReactivePower(0, 0) };
		mPower(1, 1) = { mActivePower(1, 1), mReactivePower(1, 1) };
		mPower(2, 2) = { mActivePower(2, 2), mReactivePower(2, 2) };

		mNomVoltage = std::abs(mTerminals[0]->initialSingleVoltage());

		mSLog->info("\nActive Power [W]: {}"
					"\nReactive Power [VAr]: {}",
					Logger::matrixToString(mActivePower),
					Logger::matrixToString(mReactivePower));
		mSLog->info("Nominal Voltage={} [V]", mNomVoltage);
	}

	if (mActivePower(0,0) != 0) {
		mResistance = std::pow(mNomVoltage/sqrt(3), 2) * mActivePower.inverse();
		mConductance = mResistance.inverse();
		mSubResistor = std::make_shared<EMT::Ph3::Resistor>(mName + "_res", mLogLevel);
		mSubResistor->setParameters(mResistance);
		mSubResistor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubResistor->initialize(mFrequencies);
		mSubResistor->initializeFromNodesAndTerminals(frequency);
	}

	if (mReactivePower(0, 0) != 0)
		mReactance = std::pow(mNomVoltage/sqrt(3), 2) * mReactivePower.inverse();
	else
		mReactance = Matrix::Zero(1, 1);

	if (mReactance(0,0) > 0) {
		mInductance = mReactance / (2 * PI * frequency);

		mSubInductor = std::make_shared<EMT::Ph3::Inductor>(mName + "_ind", mLogLevel);
		mSubInductor->setParameters(mInductance);
		mSubInductor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubInductor->initialize(mFrequencies);
		mSubInductor->initializeFromNodesAndTerminals(frequency);
	}
	else if (mReactance(0,0) < 0) {
		mCapacitance = -1 / (2 * PI * frequency) * mReactance.inverse();

		mSubCapacitor = std::make_shared<EMT::Ph3::Capacitor>(mName + "_cap", mLogLevel);
		mSubCapacitor->setParameters(mCapacitance);
		mSubCapacitor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubCapacitor->initialize(mFrequencies);
		mSubCapacitor->initializeFromNodesAndTerminals(frequency);
	}

	MatrixComp vInitABC = MatrixComp::Zero(3, 1);
	vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * mTerminals[0]->initialSingleVoltage();
	vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
	vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
	mIntfVoltage = vInitABC.real();

	MatrixComp iInitABC = MatrixComp::Zero(3, 1);
	// v i^T* = S
	// v^T v i^T* = v^T S
	// i^T*= (|v|^2)^(-1) v^T S

	Complex v_ = vInitABC(0, 0)*vInitABC(0, 0) + vInitABC(1, 0)*vInitABC(1, 0) + vInitABC(2, 0)*vInitABC(2, 0);
	MatrixComp rhs_ = Complex(1, 0) / v_ * vInitABC.transpose() * mPower;
	iInitABC = rhs_.conjugate().transpose();
	mIntfCurrent = iInitABC.real();

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
		Logger::matrixToString(mIntfVoltage),
		Logger::matrixToString(mIntfCurrent),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
		Logger::matrixToString(mActivePower),
		Logger::matrixToString(mReactivePower),
		Logger::matrixToString(mResistance),
		Logger::matrixToString(mReactance));
	mSLog->flush();

}

void EMT::Ph3::RXLoad::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();
	mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
	if (mSubResistor) {
		mSubResistor->mnaInitialize(omega, timeStep, leftVector);
		for (auto task : mSubResistor->mnaTasks()) {
			mMnaTasks.push_back(task);
		}
	}
	if (mSubInductor) {
		mSubInductor->mnaInitialize(omega, timeStep, leftVector);
		for (auto task : mSubInductor->mnaTasks()) {
			mMnaTasks.push_back(task);
		}
	}
	if (mSubCapacitor) {
		mSubCapacitor->mnaInitialize(omega, timeStep, leftVector);
		for (auto task : mSubCapacitor->mnaTasks()) {
			mMnaTasks.push_back(task);
		}
	}
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

void EMT::Ph3::RXLoad::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	if (mSubResistor)
		mSubResistor->mnaApplyRightSideVectorStamp(rightVector);
	if (mSubInductor)
		mSubInductor->mnaApplyRightSideVectorStamp(rightVector);
	if (mSubCapacitor)
		mSubCapacitor->mnaApplyRightSideVectorStamp(rightVector);
}

void EMT::Ph3::RXLoad::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (mSubResistor)
		mSubResistor->mnaApplySystemMatrixStamp(systemMatrix);
	if (mSubInductor)
		mSubInductor->mnaApplySystemMatrixStamp(systemMatrix);
	if (mSubCapacitor)
		mSubCapacitor->mnaApplySystemMatrixStamp(systemMatrix);
}

void EMT::Ph3::RXLoad::MnaPreStep::execute(Real time, Int timeStepCount) {
	mLoad.mnaApplyRightSideVectorStamp(mLoad.mRightVector);
}

void EMT::Ph3::RXLoad::MnaPostStep::execute(Real time, Int timeStepCount) {
	mLoad.mnaUpdateVoltage(*mLeftVector);
	mLoad.mnaUpdateCurrent(*mLeftVector);
}

void EMT::Ph3::RXLoad::mnaUpdateVoltage(const Matrix& leftVector) {
	mIntfVoltage = Matrix::Zero(3, 1);
	mIntfVoltage(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
	mIntfVoltage(1, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
	mIntfVoltage(2, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
}

void EMT::Ph3::RXLoad::mnaUpdateCurrent(const Matrix& leftVector) {
	mIntfCurrent = Matrix::Zero(3, 1);
	if (mSubResistor)
		mIntfCurrent += mSubResistor->intfCurrent();
	if (mSubInductor)
		mIntfCurrent += mSubInductor->intfCurrent();
	if (mSubCapacitor)
		mIntfCurrent += mSubCapacitor->intfCurrent();
}


// #### DAE functions ####

void EMT::Ph3::RXLoad::daeInitialize(double time, double state[], 
	double dstate_dt[], int& offset) {
	// offset: number of component in state, dstate_dt
	// state[offset] = current through voltage source flowing into node matrixNodeIndex(1)
	// dstate_dt[offset] = derivative of current through voltage source  (not used yed) --> set to zero

	updateMatrixNodeIndices();

	if (mReactance(0,0) > 0) {
		// state variable is the inductor current
		Matrix inductorCurrent = mSubInductor->intfCurrent();
		Matrix inductorVoltage = mSubInductor->intfVoltage();
		state[offset] = inductorCurrent(0,0);
		dstate_dt[offset]   = inductorVoltage(0,0)/mInductance(0,0);
		state[offset+1] = inductorCurrent(1,0);
		dstate_dt[offset+1] = inductorVoltage(1,0)/mInductance(1,1);
		state[offset+2] = inductorCurrent(2,0);
		dstate_dt[offset+2] = inductorVoltage(2,0)/mInductance(2, 2);

		mSLog->info(
			"\n--- daeInitialize ---"
			"\nmReactance(0,0) > 0  --> state variable are inductor currents"
			"\nAdded current-phase1 through the inductor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded current-phase2 through the inductor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded current-phase3 through the inductor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded derivative of current-phase1 through the inductor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\nAdded derivative of current-phase2 through the inductor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\nAdded derivative of current-phase3 through the inductor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\n--- daeInitialize finished ---",
			this->name(), state[offset],
			this->name(), state[offset+1],
			this->name(), state[offset+2],
			this->name(), dstate_dt[offset],
			this->name(), dstate_dt[offset+1],
			this->name(), dstate_dt[offset+2]
		);
	}
	else if (mReactance(0,0) < 0) {
		// state variable is the voltage through capacitor
		Matrix capacitorCurrent = mSubCapacitor->intfCurrent();
		Matrix capacitorVoltage = mSubCapacitor->intfVoltage();
		state[offset] = capacitorVoltage(0,0);
		dstate_dt[offset]   = capacitorCurrent(0,0)/mCapacitance(0,0);
		state[offset+1] = capacitorVoltage(1,0);
		dstate_dt[offset+1] = capacitorCurrent(1,0)/mCapacitance(1,1);
		state[offset+2] = capacitorVoltage(2,0);
		dstate_dt[offset+2] = capacitorCurrent(2,0)/mCapacitance(2, 2);
		
		mSLog->info(
			"\n--- daeInitialize ---"
			"\nmReactance(0,0) < 0  --> state variable are capacitor currentvoltages"
			"\nAdded voltage-phase1 through the capacitor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded voltage-phase2 through the capacitor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded voltage-phase3 through the capacitor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded derivative of voltage-phase1 through the capacitor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\nAdded derivative of voltage-phase2 through the capacitor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\nAdded derivative of voltage-phase3 through the capacitor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\n--- daeInitialize finished ---",
			this->name(), state[offset],
			this->name(), state[offset+1],
			this->name(), state[offset+2],
			this->name(), dstate_dt[offset],
			this->name(), dstate_dt[offset+1],
			this->name(), dstate_dt[offset+2]
		);
	}
	mSLog->flush();
	offset+=3;
}

void EMT::Ph3::RXLoad::daeResidual(double sim_time, 
	const double state[], const double dstate_dt[], 
	double resid[], std::vector<int>& off) {
	/*
	v = R*i + L*der(i);
	*/

	int c_offset = off[0]+off[1]; //current offset for component

	if (mReactance(0,0) > 0) {
		resid[c_offset]   = state[matrixNodeIndex(0, 0)] - mInductance(0,0)*dstate_dt[c_offset]   - state[c_offset]*mResistance(0,0);
		resid[c_offset+1] = state[matrixNodeIndex(0, 1)] - mInductance(1,1)*dstate_dt[c_offset+1] - state[c_offset+1]*mResistance(1,1);
		resid[c_offset+2] = state[matrixNodeIndex(0, 2)] - mInductance(2,2)*dstate_dt[c_offset+2] - state[c_offset+2]*mResistance(2,2);
		resid[matrixNodeIndex(0, 0)] += state[c_offset] ;
		resid[matrixNodeIndex(0, 1)] += state[c_offset+1];
		resid[matrixNodeIndex(0, 2)] += state[c_offset+2];

	}
	else if (mReactance(0,0) < 0)
	{
		resid[c_offset]   = (state[matrixNodeIndex(0, 0)]-state[c_offset])*mConductance(0,0)   - mCapacitance(0,0)*dstate_dt[c_offset];
		resid[c_offset+1] = (state[matrixNodeIndex(0, 1)]-state[c_offset+1])*mConductance(1,1) - mCapacitance(0,0)*dstate_dt[c_offset+1];
		resid[c_offset+2] = (state[matrixNodeIndex(0, 2)]-state[c_offset+2])*mConductance(2,2) - mCapacitance(0,0)*dstate_dt[c_offset+2];
		resid[matrixNodeIndex(0, 0)] += (state[matrixNodeIndex(0, 0)] - state[c_offset])*mConductance(0,0);
		resid[matrixNodeIndex(0, 1)] += (state[matrixNodeIndex(0, 1)] - state[c_offset+1])*mConductance(1,1);
		resid[matrixNodeIndex(0, 2)] += (state[matrixNodeIndex(0, 2)] - state[c_offset+2])*mConductance(2,2);
	}
	
	off[1] += 3;
}

void EMT::Ph3::RXLoad::daePostStep(const double state[], const double dstate_dt[], int& offset) {
	mIntfVoltage(0, 0) = state[matrixNodeIndex(0, 0)];
	mIntfVoltage(1, 0) = state[matrixNodeIndex(0, 1)];
	mIntfVoltage(2, 0) = state[matrixNodeIndex(0, 2)];
	if (mReactance(0,0) > 0) {
		mIntfCurrent(0, 0) = state[offset++];
		mIntfCurrent(1, 0) = state[offset++];
		mIntfCurrent(2, 0) = state[offset++];
	}
	else if (mReactance(0,0) < 0)
	{
		mIntfCurrent(0, 0) = (state[matrixNodeIndex(0, 0)]-state[offset++])*mConductance(0,0);
		mIntfCurrent(1, 0) = (state[matrixNodeIndex(0, 1)]-state[offset++])*mConductance(1,1);
		mIntfCurrent(2, 0) = (state[matrixNodeIndex(0, 2)]-state[offset++])*mConductance(2,2);
	}
}

