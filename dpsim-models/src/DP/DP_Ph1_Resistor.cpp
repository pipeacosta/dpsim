/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/DP/DP_Ph1_Resistor.h>

using namespace CPS;

DP::Ph1::Resistor::Resistor(String uid, String name, Logger::Level logLevel)
	: Base::Ph1::Resistor(mAttributes), SimPowerComp<Complex>(uid, name, logLevel) {
	**mIntfVoltage = MatrixComp::Zero(1,1);
	**mIntfCurrent = MatrixComp::Zero(1,1);
	setTerminalNumber(2);
}

SimPowerComp<Complex>::Ptr DP::Ph1::Resistor::clone(String name) {
	auto copy = Resistor::make(name, mLogLevel);
	copy->setParameters(**mResistance);
	return copy;
}

void DP::Ph1::Resistor::initializeFromNodesAndTerminals(Real frequency) {

	Complex impedance = { **mResistance, 0 };
	(**mIntfVoltage)(0,0) = initialSingleVoltage(1) - initialSingleVoltage(0);
	(**mIntfCurrent)(0,0) = (**mIntfVoltage)(0,0) / impedance;

	mSLog->info("\nResistance [Ohm]: {:s}"
				"\nImpedance [Ohm]: {:s}",
				Logger::realToString(**mResistance),
				Logger::complexToString(impedance));
	mSLog->info("\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nTerminal 1 voltage: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::phasorToString((**mIntfVoltage)(0,0)),
		Logger::phasorToString((**mIntfCurrent)(0,0)),
		Logger::phasorToString(initialSingleVoltage(0)),
		Logger::phasorToString(initialSingleVoltage(1)));
}

// #### MNA functions ####
void DP::Ph1::Resistor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	mSLog->info(
		"\n--- MNA initialization ---"
		"\nInitial voltage {:s}"
		"\nInitial current {:s}"
		"\n--- MNA initialization finished ---",
		Logger::phasorToString((**mIntfVoltage)(0, 0)),
		Logger::phasorToString((**mIntfCurrent)(0, 0)));
}

void DP::Ph1::Resistor::mnaInitializeHarm(Real omega, Real timeStep, std::vector<Attribute<Matrix>::Ptr> leftVectors) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStepHarm>(*this, leftVectors));
}

void DP::Ph1::Resistor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	Complex conductance = Complex(1. / **mResistance, 0);

	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		// Set diagonal entries
		if (terminalNotGrounded(0))
			Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), conductance, mNumFreqs, freq);
		if (terminalNotGrounded(1))
		// Set off diagonal entries
			Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), conductance, mNumFreqs, freq);
		if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
			Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -conductance, mNumFreqs, freq);
			Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -conductance, mNumFreqs, freq);
		}

		mSLog->info("-- Stamp frequency {:d} ---", freq);
		if (terminalNotGrounded(0))
			mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(conductance), matrixNodeIndex(0), matrixNodeIndex(0));
		if (terminalNotGrounded(1))
			mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(conductance), matrixNodeIndex(1), matrixNodeIndex(1));
		if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
			mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(-conductance), matrixNodeIndex(0), matrixNodeIndex(1));
			mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(-conductance), matrixNodeIndex(1), matrixNodeIndex(0));
		}
	}
}

void DP::Ph1::Resistor::mnaApplySystemMatrixStampHarm(Matrix& systemMatrix, Int freqIdx) {
	Complex conductance = Complex(1. / **mResistance, 0);
	// Set diagonal entries
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), conductance);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), conductance);
	// Set off diagonal entries
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -conductance);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -conductance);
	}

	mSLog->info("-- Stamp for frequency {:f} ---", mFrequencies(freqIdx,0));
	if (terminalNotGrounded(0))
		mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(conductance), matrixNodeIndex(0), matrixNodeIndex(0));
	if (terminalNotGrounded(1))
		mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(conductance), matrixNodeIndex(1), matrixNodeIndex(1));
	if (terminalNotGrounded(0)  &&  terminalNotGrounded(1)) {
		mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(-conductance), matrixNodeIndex(0), matrixNodeIndex(1));
		mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(-conductance), matrixNodeIndex(1), matrixNodeIndex(0));
	}
}

void DP::Ph1::Resistor::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(this->mIntfVoltage);
	modifiedAttributes.push_back(this->mIntfCurrent);
}

void DP::Ph1::Resistor::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	this->mnaUpdateVoltage(**leftVector);
	this->mnaUpdateCurrent(**leftVector);
}

void DP::Ph1::Resistor::MnaPostStepHarm::execute(Real time, Int timeStepCount) {
	for (UInt freq = 0; freq < mResistor.mNumFreqs; freq++)
		mResistor.mnaUpdateVoltageHarm(**mLeftVectors[freq], freq);
	mResistor.mnaUpdateCurrentHarm();
}

void DP::Ph1::Resistor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		(**mIntfVoltage)(0,freq) = 0;
		if (terminalNotGrounded(1))
			(**mIntfVoltage)(0,freq) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1), mNumFreqs, freq);
		if (terminalNotGrounded(0))
			(**mIntfVoltage)(0,freq) = (**mIntfVoltage)(0,freq) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0), mNumFreqs, freq);

		SPDLOG_LOGGER_DEBUG(mSLog, "Voltage {:s}", Logger::phasorToString((**mIntfVoltage)(0,freq)));
	}
}

void DP::Ph1::Resistor::mnaUpdateCurrent(const Matrix& leftVector) {
	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		(**mIntfCurrent)(0,freq) = (**mIntfVoltage)(0,freq) / **mResistance;
		SPDLOG_LOGGER_DEBUG(mSLog, "Current {:s}", Logger::phasorToString((**mIntfCurrent)(0,freq)));
	}
}

void DP::Ph1::Resistor::mnaUpdateVoltageHarm(const Matrix& leftVector, Int freqIdx) {
	// v1 - v0
	(**mIntfVoltage)(0,freqIdx) = 0;
	if (terminalNotGrounded(1))
		(**mIntfVoltage)(0,freqIdx) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1));
	if (terminalNotGrounded(0))
		(**mIntfVoltage)(0,freqIdx) = (**mIntfVoltage)(0,freqIdx) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0));

	SPDLOG_LOGGER_DEBUG(mSLog, "Voltage {:s}", Logger::phasorToString((**mIntfVoltage)(0,freqIdx)));
}

void DP::Ph1::Resistor::mnaUpdateCurrentHarm() {
	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		(**mIntfCurrent)(0,freq) = (**mIntfVoltage)(0,freq) / **mResistance;
		SPDLOG_LOGGER_DEBUG(mSLog, "Current {:s}", Logger::phasorToString((**mIntfCurrent)(0,freq)));
	}
}

void DP::Ph1::Resistor::mnaTearApplyMatrixStamp(Matrix& tearMatrix) {
	Math::addToMatrixElement(tearMatrix, mTearIdx, mTearIdx, Complex(**mResistance, 0));
}

// #### DAE functions ####

void DP::Ph1::Resistor::daeResidual(double ttime, const double state[], const double dstate_dt[], double resid[], std::vector<int>& off) {
	// new state vector definintion:
	// state[0]=node0_voltage
	// state[1]=node1_voltage
	// ....
	// state[n]=noden_voltage
	// state[n+1]=component0_voltage
	// state[n+2]=component0_inductance (not yet implemented)
	// ...
	// state[m-1]=componentm_voltage
	// state[m]=componentm_inductance
	// ...
	// state[x] = nodal_equation_1
	// state[x+1] = nodal_equation_2
	// ...

    int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);
	int c_offset = off[0]+off[1]; //current offset for component
	int n_offset_1 = c_offset + Pos1 + 1;// current offset for first nodal equation
	int n_offset_2 = c_offset + Pos2 + 1;// current offset for second nodal equation
	resid[c_offset] = (state[Pos2]-state[Pos1]) - state[c_offset]; // Voltage equation for Resistor
	//resid[c_offset+1] = ; //TODO : add inductance equation
    resid[n_offset_1] += 1.0/ **mResistance * state[c_offset];
    resid[n_offset_2] += 1.0/ **mResistance * state[c_offset];
	off[1] += 1;
}

Complex DP::Ph1::Resistor::daeInitialize() {

	 return (**mIntfVoltage)(0,0);
}
