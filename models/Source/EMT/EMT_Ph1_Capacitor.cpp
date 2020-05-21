/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/EMT/EMT_Ph1_Capacitor.h>

using namespace CPS;

EMT::Ph1::Capacitor::Capacitor(String uid, String name,	Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	mEquivCurrent = 0;
	mIntfVoltage = Matrix::Zero(1,1);
	mIntfCurrent = Matrix::Zero(1,1);
	setTerminalNumber(2);

	addAttribute<Real>("C", &mCapacitance, Flags::read | Flags::write);
}

SimPowerComp<Real>::Ptr EMT::Ph1::Capacitor::clone(String name) {
	auto copy = Capacitor::make(name, mLogLevel);
	copy->setParameters(mCapacitance);
	return copy;
}

void EMT::Ph1::Capacitor::initializeFromNodesAndTerminals(Real frequency) {

	Real omega = 2 * PI * frequency;
	Complex impedance = { 0, - 1. / (omega * mCapacitance) };
	mIntfVoltage(0,0) = (initialSingleVoltage(1) - initialSingleVoltage(0)).real();
	mIntfCurrent(0,0) = ((initialSingleVoltage(1) - initialSingleVoltage(0)) / impedance).real();

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

void EMT::Ph1::Capacitor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mEquivCond = (2.0 * mCapacitance) / timeStep;
	// Update internal state
	mEquivCurrent = -mIntfCurrent(0,0) + -mEquivCond * mIntfVoltage(0,0);

	mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

void EMT::Ph1::Capacitor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), mEquivCond);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), mEquivCond);
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -mEquivCond);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -mEquivCond);
	}
}

void EMT::Ph1::Capacitor::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	mEquivCurrent = -mIntfCurrent(0,0) + -mEquivCond * mIntfVoltage(0,0);
	if (terminalNotGrounded(0))
		Math::setVectorElement(rightVector, matrixNodeIndex(0), mEquivCurrent);
	if (terminalNotGrounded(1))
		Math::setVectorElement(rightVector, matrixNodeIndex(1), -mEquivCurrent);
}

void EMT::Ph1::Capacitor::MnaPreStep::execute(Real time, Int timeStepCount) {
	mCapacitor.mnaApplyRightSideVectorStamp(mCapacitor.mRightVector);
}

void EMT::Ph1::Capacitor::MnaPostStep::execute(Real time, Int timeStepCount) {
	mCapacitor.mnaUpdateVoltage(*mLeftVector);
	mCapacitor.mnaUpdateCurrent(*mLeftVector);
}

void EMT::Ph1::Capacitor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	mIntfVoltage(0,0) = 0;
	if (terminalNotGrounded(1))
		mIntfVoltage(0,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1));
	if (terminalNotGrounded(0))
		mIntfVoltage(0,0) = mIntfVoltage(0,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0));
}

void EMT::Ph1::Capacitor::mnaUpdateCurrent(const Matrix& leftVector) {
	mIntfCurrent(0,0) = mEquivCond * mIntfVoltage(0,0) + mEquivCurrent;
}


// #### DAE functions ####


void EMT::Ph1::Capacitor::daeInitialize(double time, double state[], double dstate_dt[], int &offset) {
	// state[offset] = voltage through capacitor
	// dstate_dt[offset] = voltage derivative though capacitor 
	updateMatrixNodeIndices();
	state[offset] = mIntfVoltage(0,0);
	dstate_dt[offset] = mIntfCurrent(0,0)/mCapacitance;		

	mSLog->info(
		"\n--- daeInitialize ---",
		"\nAdded voltage through capacitor '{:s}' to state vector, initial value  = {:f}V", 
		"\nAdded derivative of voltage through capacitor '{:s}' to derivative state vector, initial value= {:f}",
		"\n--- daeInitialize finished ---",
		this->name(),state[offset],
		this->name(),dstate_dt[offset]
	);

	offset++;
}

void EMT::Ph1::Capacitor::daeResidual(double sim_time, 
	const double state[], const double dstate_dt[], 
	double resid[], std::vector<int>& off) {
	// state[c_offset] = voltage through capacitor
	// dstate_dt[c_offset] = voltage capacitor derivative
	// state[Pos2] = voltage of node matrixNodeIndex(1)
	// state[Pos1] = voltage of node matrixNodeIndex(0)
	// resid[c_offset] = voltage eq of capacitor: state[pos2]-state[pos1]-state[offset]=0
	// resid[Pos1] = nodal current equation of node matrixNodeIndex(0) ---> substract current through cap
	// resid[Pos2] = nodal current equation of node matrixNodeIndex(1) ---> add current through cap

	int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);
	int c_offset = off[0]+off[1]; //current offset for component

	resid[c_offset] = -state[c_offset];
	if (terminalNotGrounded(0)) {
		resid[Pos1] -= mCapacitance*dstate_dt[c_offset];
		resid[c_offset] -= state[Pos1];
	}
	if (terminalNotGrounded(1)) {
		resid[Pos2] += mCapacitance*dstate_dt[c_offset];
		resid[c_offset] += state[Pos2];
	}

	off[1] += 1;
}

void EMT::Ph1::Capacitor::daePostStep(const double state[], const double dstate_dt[], int& offset) {
	mIntfCurrent(0,0) = mCapacitance*dstate_dt[offset];
	mIntfVoltage(0,0) = state[offset];
	offset++;
}

