/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/EMT/EMT_Ph1_Resistor.h>

using namespace CPS;

EMT::Ph1::Resistor::Resistor(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	setTerminalNumber(2);
	mIntfVoltage = Matrix::Zero(1,1);
	mIntfCurrent = Matrix::Zero(1,1);

	addAttribute<Real>("R", &mResistance, Flags::read | Flags::write);
}

SimPowerComp<Real>::Ptr EMT::Ph1::Resistor::clone(String name) {
	auto copy = Resistor::make(name, mLogLevel);
	copy->setParameters(mResistance);
	return copy;
}

void EMT::Ph1::Resistor::initializeFromNodesAndTerminals(Real frequency) {

	mIntfVoltage(0,0) = (initialSingleVoltage(1) - initialSingleVoltage(0)).real();
	mIntfCurrent(0,0) = mIntfVoltage(0,0) / mResistance;

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

// #### MNA functions ####
void EMT::Ph1::Resistor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

void EMT::Ph1::Resistor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	Real conductance = 1. / mResistance;
	// Set diagonal entries
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), conductance);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), conductance);
	// Set off diagonal entries
	if (terminalNotGrounded(0)  &&  terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -conductance);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -conductance);
	}

	if (terminalNotGrounded(0))
		mSLog->info("Add {:f} to system at ({:d},{:d})", conductance, matrixNodeIndex(0), matrixNodeIndex(0));
	if (terminalNotGrounded(1))
		mSLog->info("Add {:f} to system at ({:d},{:d})", conductance, matrixNodeIndex(1), matrixNodeIndex(1));
	if ( terminalNotGrounded(0)  &&  terminalNotGrounded(1) ) {
		mSLog->info("Add {:f} to system at ({:d},{:d})", -conductance, matrixNodeIndex(0), matrixNodeIndex(1));
		mSLog->info("Add {:f} to system at ({:d},{:d})", -conductance, matrixNodeIndex(1), matrixNodeIndex(0));
	}
}

void EMT::Ph1::Resistor::MnaPostStep::execute(Real time, Int timeStepCount) {
	mResistor.mnaUpdateVoltage(*mLeftVector);
	mResistor.mnaUpdateCurrent(*mLeftVector);
}

void EMT::Ph1::Resistor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	mIntfVoltage(0,0) = 0;
	if (terminalNotGrounded(1))
		mIntfVoltage(0,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1));
	if (terminalNotGrounded(0))
		mIntfVoltage(0,0) = mIntfVoltage(0,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0));
}

void EMT::Ph1::Resistor::mnaUpdateCurrent(const Matrix& leftVector) {
	mIntfCurrent(0,0) = mIntfVoltage(0,0) / mResistance;
	mSLog->info(
		"\n--- mnaUpdateCurrent ---"
		"\nmIntfCurrent(0,0): {:f}"
		"\n--- mnaUpdateCurrent finished ---",
		mIntfCurrent(0,0));
}

// #### DAE functions ####

void EMT::Ph1::Resistor::daeInitialize(double time, double state[], double dstate_dt[], int& offset) {
	updateMatrixNodeIndices();
	mSLog->info(
		"\n--- daeInitialize ---");
}

void EMT::Ph1::Resistor::daeResidual(double sim_time,
	const double state[], const double dstate_dt[],
	double resid[], std::vector<int>& off) {
	// state[Pos2] = voltage of node matrixNodeIndex(1)
	// state[Pos1] = voltage of node matrixNodeIndex(0)
	// resid[Pos1] = nodal current equation of node matrixNodeIndex(1) --> add resistor current
	// resid[Pos2] = nodal current equation of node matrixNodeIndex(0) --> substract resistor current

	int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);

	double voltageResistor = 0.0;
	if (terminalNotGrounded(0)) {
		voltageResistor -= state[Pos1];
	}
	if (terminalNotGrounded(1)) {
		voltageResistor +=state[Pos2];
	}
	if (terminalNotGrounded(0)) {
		resid[Pos1] -= voltageResistor/mResistance;
	}
	if (terminalNotGrounded(1)) {
		resid[Pos2] += voltageResistor/mResistance;
	}
}

void EMT::Ph1::Resistor::daePostStep(const double state[], const double dstate_dt[], int& offset) {
	int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);
	mIntfVoltage(0,0) = 0.0;
	if (terminalNotGrounded(0)) {
		mIntfVoltage(0,0) -= state[Pos1];
	}
	if (terminalNotGrounded(1)) {
		mIntfVoltage(0,0) += state[Pos2];
	}
	mIntfCurrent(0,0) = mIntfVoltage(0,0) / mResistance;
}
