/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/EMT/EMT_Ph1_Inductor.h>

using namespace CPS;

EMT::Ph1::Inductor::Inductor(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	mEquivCurrent = 0;
	mIntfVoltage = Matrix::Zero(1,1);
	mIntfCurrent = Matrix::Zero(1,1);
	setTerminalNumber(2);

	addAttribute<Real>("L", &mInductance, Flags::read | Flags::write);
}

SimPowerComp<Real>::Ptr EMT::Ph1::Inductor::clone(String name) {
	auto copy = Inductor::make(name, mLogLevel);
	copy->setParameters(mInductance);
	return copy;
}

void EMT::Ph1::Inductor::initializeFromNodesAndTerminals(Real frequency) {

	Real omega = 2 * PI * frequency;
	Complex impedance = { 0, omega * mInductance };
	mIntfVoltage(0,0) = (initialSingleVoltage(1) - initialSingleVoltage(0)).real();
	mIntfCurrent(0,0) = (mIntfVoltage(0,0) / impedance).real();

	mSLog->info(
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:f}"
		"\nCurrent: {:f}"
		"\nTerminal 0 voltage: {:f}"
		"\nTerminal 1 voltage: {:f}"
		"\n--- Initialization from powerflow finished ---",
		mIntfVoltage(0,0),
		mIntfCurrent(0,0),
		initialSingleVoltage(0).real(),
		initialSingleVoltage(1).real());
}

void EMT::Ph1::Inductor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mEquivCond = timeStep / (2.0 * mInductance);
	// Update internal state
	mEquivCurrent = mEquivCond * mIntfVoltage(0,0) + mIntfCurrent(0,0);
	mSLog->info(
		"\n--- mnaInitialize ---"
		"\nmEquivCond: {:f}"
		"\nmIntfVoltage(0,0): {:f}"
		"\nmIntfCurrent(0,0): {:f}"
		"\nmEquivCurrent: {:f}"
		"\n--- Initialization mnaInitialize finished ---",
		mEquivCond,
		mIntfCurrent(0,0),
		mIntfVoltage(0,0),
		mIntfCurrent(0,0),
		mEquivCurrent);

	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void EMT::Ph1::Inductor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), mEquivCond);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), mEquivCond);
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -mEquivCond);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -mEquivCond);
	}
}

void EMT::Ph1::Inductor::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	// Update internal state
	mEquivCurrent = mEquivCond * mIntfVoltage(0,0) + mIntfCurrent(0,0);
	if (terminalNotGrounded(0))
		Math::setVectorElement(rightVector, matrixNodeIndex(0), mEquivCurrent);
	if (terminalNotGrounded(1))
		Math::setVectorElement(rightVector, matrixNodeIndex(1), -mEquivCurrent);
}

void EMT::Ph1::Inductor::MnaPreStep::execute(Real time, Int timeStepCount) {
	mInductor.mnaApplyRightSideVectorStamp(mInductor.mRightVector);
}

void EMT::Ph1::Inductor::MnaPostStep::execute(Real time, Int timeStepCount) {
	mInductor.mnaUpdateVoltage(*mLeftVector);
	mInductor.mnaUpdateCurrent(*mLeftVector);
}

void EMT::Ph1::Inductor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	mIntfVoltage(0,0) = 0;
	if (terminalNotGrounded(1))
		mIntfVoltage(0,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1));
	if (terminalNotGrounded(0))
		mIntfVoltage(0,0) = mIntfVoltage(0,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0));
}

void EMT::Ph1::Inductor::mnaUpdateCurrent(const Matrix& leftVector) {
	mIntfCurrent(0,0) = mEquivCond * mIntfVoltage(0,0) + mEquivCurrent;
	mSLog->info(
		"\n--- mnaUpdateCurrent ---"
		"\nmEquivCond: {:f}"
		"\nmIntfVoltage(0,0): {:f}"
		"\nmIntfCurrent(0,0): {:f}"
		"\nmEquivCurrent: {:f}"
		"\n--- Initialization mnaUpdateCurrent finished ---",
		mEquivCond,
		mIntfCurrent(0,0),
		mIntfVoltage(0,0),
		mIntfCurrent(0,0),
		mEquivCurrent);
}

// #### DAE functions ####


void EMT::Ph1::Inductor::daeInitialize(double state[], double dstate_dt[], int &counter) {
	// state[c_offset] = current through inductor
	// dstate_dt[c_offset] = inductor current derivative

	updateMatrixNodeIndices();
	
	int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);
	mIntfVoltage(0,0) = 0.0;
	if (terminalNotGrounded(0)) {
		mIntfVoltage(0,0) -= state[Pos1];
	}
	if (terminalNotGrounded(1)) {
		mIntfVoltage(0,0) += state[Pos2];
	}
	mIntfCurrent(0,0) = 0.0;											//TODO

	state[counter] = mIntfCurrent(0,0);	
	dstate_dt[counter] = mIntfVoltage(0,0)/mInductance;					//correct?
	mSLog->info("mIntfVoltage(0,0) = {}V", mIntfVoltage(0,0));
	mSLog->info("Added current through inductor '{:s}' to state vector, initial value = {}A", this->name(), state[counter]);
	mSLog->info("Added derivative of current inductor '{:s}' to derivative state vector, initial value = {}A", this->name(), dstate_dt[counter]);
	counter++;
}

void EMT::Ph1::Inductor::daeResidual(double ttime, const double state[], const double dstate_dt[], double resid[], std::vector<int>& off) {
	// state[c_offset] = current through inductor, flowing into node matrixNodeIndex(0)
	// dstate_dt[c_offset] = inductor current derivative
	// state[Pos2] = voltage of node matrixNodeIndex(1)
	// state[Pos1] = voltage of node matrixNodeIndex(0)
	// resid[c_offset] = voltage eq of inductor: v_ind(t) -L*(d/dt)i_ind(t) = 0 --> state[Pos2] - state[Pos1] - L*dstate_dt[c_offset+1] = 0
	// resid[Pos1] = nodal current equation of node matrixNodeIndex(0), flowing into node matrixNodeIndex(0) ---> substract state[c_offset+1]
	// resid[Pos2] = nodal current equation of node matrixNodeIndex(1), flowing into node matrixNodeIndex(0) ---> add state[c_offset+1]

	int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);
	int c_offset = off[0]+off[1]; //current offset for component

	resid[c_offset] = -mInductance*dstate_dt[c_offset];
	if (terminalNotGrounded(0)) {
		resid[c_offset] -= state[Pos1];
		resid[Pos1] -= state[c_offset];
	}
	if (terminalNotGrounded(1)) {
		resid[c_offset] += state[Pos2];
		resid[Pos2] += state[c_offset];
	}

	mSLog->info("state[{}]={}, Pos2={}", Pos2, state[Pos2], Pos2);
	mSLog->info("state[c_offset]={}", state[c_offset]);
	mSLog->info("dstate[c_offset]={}", dstate_dt[c_offset]);
	mSLog->info("resid[{}]={}", c_offset, resid[c_offset]);
	mSLog->info("");
	off[1] += 1;
}

void EMT::Ph1::Inductor::daePostStep(const double state[], int& counter, double time) {
	int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);
	mIntfVoltage(0,0) = 0.0;
	if (terminalNotGrounded(0)) {
		mIntfVoltage(0,0) -= state[Pos1];
	}
	if (terminalNotGrounded(1)) {
		mIntfVoltage(0,0) += state[Pos2];
	}
	mIntfCurrent(0,0) = state[counter];
	counter++;
}