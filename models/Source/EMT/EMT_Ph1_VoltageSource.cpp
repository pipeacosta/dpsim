/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/EMT/EMT_Ph1_VoltageSource.h>

using namespace CPS;

EMT::Ph1::VoltageSource::VoltageSource(String uid, String name,	Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	setVirtualNodeNumber(1);
	setTerminalNumber(2);
	mIntfVoltage = Matrix::Zero(1,1);
	mIntfCurrent = Matrix::Zero(1,1);

	addAttribute<Complex>("V_ref", Flags::read | Flags::write);
	addAttribute<Real>("f_src", Flags::read | Flags::write);
}

void EMT::Ph1::VoltageSource::setParameters(Complex voltageRef, Real srcFreq) {
	attribute<Complex>("V_ref")->set(voltageRef);
	attribute<Real>("f_src")->set(srcFreq);

	mParametersSet = true;
}

SimPowerComp<Real>::Ptr EMT::Ph1::VoltageSource::clone(String name) {
	auto copy = VoltageSource::make(name, mLogLevel);
	copy->setParameters(attribute<Complex>("V_ref")->get(), attribute<Real>("f_src")->get());
	return copy;
}

void EMT::Ph1::VoltageSource::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mVoltageRef = attribute<Complex>("V_ref");
	mSrcFreq = attribute<Real>("f_src");
	mIntfVoltage(0,0) = Math::abs(mVoltageRef->get()) * cos(Math::phase(mVoltageRef->get()));
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void EMT::Ph1::VoltageSource::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (terminalNotGrounded(0)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex(), -1);
		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0), -1);
	}
	if (terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex(), 1);
		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1), 1);
	}

	if (terminalNotGrounded(0)) {
		mSLog->info("Add {:f} to system at ({:d},{:d})", -1, matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex());
		mSLog->info("Add {:f} to system at ({:d},{:d})", -1, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0));
	}
	if (terminalNotGrounded(1)) {
		mSLog->info("Add {:f} to system at ({:d},{:d})", 1, matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex());
		mSLog->info("Add {:f} to system at ({:d},{:d})", 1, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1));
	}
}

void EMT::Ph1::VoltageSource::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	Math::setVectorElement(rightVector, mVirtualNodes[0]->matrixNodeIndex(), mIntfVoltage(0,0));
}

void EMT::Ph1::VoltageSource::updateVoltage(Real time) {
	Complex voltageRef = mVoltageRef->get();
	Real srcFreq = mSrcFreq->get();
	if (srcFreq > 0)
		mIntfVoltage(0,0) = Math::abs(voltageRef) * cos(time * 2.*PI*srcFreq + Math::phase(voltageRef));
	else
		mIntfVoltage(0,0) = voltageRef.real();
}

void EMT::Ph1::VoltageSource::MnaPreStep::execute(Real time, Int timeStepCount) {
	mVoltageSource.updateVoltage(time);
	mVoltageSource.mnaApplyRightSideVectorStamp(mVoltageSource.mRightVector);
}

void EMT::Ph1::VoltageSource::MnaPostStep::execute(Real time, Int timeStepCount) {
	mVoltageSource.mnaUpdateCurrent(*mLeftVector);
}

void EMT::Ph1::VoltageSource::mnaUpdateCurrent(const Matrix& leftVector) {
	mIntfCurrent(0,0) = Math::realFromVectorElement(leftVector, mVirtualNodes[0]->matrixNodeIndex());
}

// #### DAE functions ####

void EMT::Ph1::VoltageSource::daeUpdateVoltage(Real time) {
	mVoltageRef = attribute<Complex>("V_ref");
	mSrcFreq = attribute<Real>("f_src");
	Complex voltageRef = mVoltageRef->get();
	Real srcFreq = mSrcFreq->get();
	if (srcFreq > 0)
		mIntfVoltage(0,0) = Math::abs(voltageRef) * cos(time * 2.*PI*srcFreq + Math::phase(voltageRef));
	else
		mIntfVoltage(0,0) = voltageRef.real();
}

void EMT::Ph1::VoltageSource::daeInitialize(double state[], double dstate_dt[], int& counter) {
	// state[c_offset] = current through voltage source flowing into node matrixNodeIndex(1)

	//update initial voltage source voltage and current
	updateMatrixNodeIndices();
	double ttime = 0.0; 				//TODO: ttime=initialTime
	this->daeUpdateVoltage(ttime);
	state[counter] = mIntfCurrent(0,0);
	dstate_dt[counter] = 0;
	mSLog->info("Initial voltage of VoltageSource '{:s}' = {}V", this->name(), mIntfVoltage(0,0));
	mSLog->info("Added current of VoltageSource '{:s}' to state vector, initial value i_{:s} = {}A", this->name(), this->name(), state[counter]);
	counter++;
}

void EMT::Ph1::VoltageSource::daeResidual(double ttime, const double state[], const double dstate_dt[], double resid[], std::vector<int>& off){
	// state[c_offset] = current through voltage source flowing into node matrixNodeIndex(1)
	// resid[c_offset] = v2-v1-v_s = state[Pos2] - state[Pos1] - mIntfVoltage(0,0)
	// resid[Pos1] = nodal current equation of node matrixNodeIndex(0)
	// resid[Pos2] = nodal current equation of node matrixNodeIndex(1)

	int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);
	int c_offset = off[0]+off[1]; //current offset for component

	//reset residual vector
	//TODO: where must be reseted?????
	resid[0] = 0.0;
	resid[1] = 0.0;
	resid[c_offset] = -mIntfVoltage(0,0);
	if (terminalNotGrounded(0)) {
		resid[c_offset] -= state[Pos1];
		resid[Pos1] += state[c_offset];
	}
	if (terminalNotGrounded(1)) {
		resid[c_offset] += state[Pos2];
		resid[Pos2] -= state[c_offset];
	}
	mSLog->info("state[Pos1]={}, Pos1={}", state[Pos1], Pos1);
	mSLog->info("state[Pos2]={}, Pos2={}", state[Pos2], Pos2);
	mSLog->info("mIntfVoltage(0,0)={}", mIntfVoltage(0,0));
	mSLog->info("resid[{}]={}", c_offset, resid[c_offset]);
	mSLog->info("");
	off[1] += 1;
}

void EMT::Ph1::VoltageSource::daePostStep(const double state[], int& counter, double time) {
	this->daeUpdateVoltage(time);
	mIntfCurrent(0,0) = state[counter];
	counter++;
}