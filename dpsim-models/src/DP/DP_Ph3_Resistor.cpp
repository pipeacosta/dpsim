/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/DP/DP_Ph3_Resistor.h>

using namespace CPS;

DP::Ph3::Resistor::Resistor(String uid, String name,
	Logger::Level logLevel)
	: Base::Ph3::Resistor(mAttributes), SimPowerComp<Complex>(uid, name, logLevel) {
	mPhaseType = PhaseType::ABC;
	setTerminalNumber(2);
	**mIntfVoltage = MatrixComp::Zero(3,1);
	**mIntfCurrent = MatrixComp::Zero(3,1);
}

SimPowerComp<Complex>::Ptr DP::Ph3::Resistor::clone(String name) {
	auto copy = Resistor::make(name, mLogLevel);
	copy->setParameters(**mResistance);
	return copy;
}

void DP::Ph3::Resistor::initializeFromNodesAndTerminals(Real frequency) {

	Real voltMag = Math::abs((**mIntfVoltage)(0, 0));
	Real voltPhase = Math::phase((**mIntfVoltage)(0, 0));
	(**mIntfVoltage)(1, 0) = Complex(
		voltMag * cos(voltPhase - 2. / 3. * M_PI),
		voltMag * sin(voltPhase - 2. / 3. * M_PI));
	(**mIntfVoltage)(2, 0) = Complex(
		voltMag * cos(voltPhase + 2. / 3. * M_PI),
		voltMag * sin(voltPhase + 2. / 3. * M_PI));
	**mIntfCurrent = (**mResistance).inverse() * **mIntfVoltage;

	mSLog->info("\n--- Initialization from powerflow ---"
		"\nVoltage across amplitude and phase: \n{:s}"
		"\nCurrent amplitude and phase: \n{:s}"
		"\nTerminal 0 voltage amplitude and phase: \n{:s}"
		"\nTerminal 1 voltage amplitude and phase: \n{:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::phasorMatrixToString(**mIntfVoltage),
		Logger::phasorMatrixToString(**mIntfCurrent),
		Logger::phasorMatrixToString(initialVoltage(0)),
		Logger::phasorMatrixToString(initialVoltage(1)));
}

void DP::Ph3::Resistor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

void DP::Ph3::Resistor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {

	Matrix conductance = (**mResistance).inverse();

	//// Set diagonal entries
	//if (terminalNotGrounded(0))
	//	Math::addToMatrixElement(systemMatrix, matrixNodeIndices(0), matrixNodeIndices(0), conductance);
	//if (terminalNotGrounded(1))
	//	Math::addToMatrixElement(systemMatrix, matrixNodeIndices(1), matrixNodeIndices(1), conductance);
	//// Set off diagonal entries
	//if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
	//	Math::addToMatrixElement(systemMatrix, matrixNodeIndices(0), matrixNodeIndices(1), -conductance);
	//	Math::addToMatrixElement(systemMatrix, matrixNodeIndices(1), matrixNodeIndices(0), -conductance);
	//}
	// Set diagonal entries
	if (terminalNotGrounded(0)) {
		// set upper left block, 3x3 entries
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0), Complex(conductance(0, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 1), Complex(conductance(0, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 2), Complex(conductance(0, 2), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 0), Complex(conductance(1, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 1), Complex(conductance(1, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 2), Complex(conductance(1, 2), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 0), Complex(conductance(2, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 1), Complex(conductance(2, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 2), Complex(conductance(2, 2), 0));
	}
	if (terminalNotGrounded(1)) {
		// set buttom right block, 3x3 entries
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 0), Complex(conductance(0, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 1), Complex(conductance(0, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 2), Complex(conductance(0, 2), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 0), Complex(conductance(1, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 1), Complex(conductance(1, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 2), Complex(conductance(1, 2), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 0), Complex(conductance(2, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 1), Complex(conductance(2, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 2), Complex(conductance(2, 2), 0));
	}
	// Set off diagonal blocks, 2x3x3 entries
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 0), -Complex(conductance(0, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 1), -Complex(conductance(0, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 2), -Complex(conductance(0, 2), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 0), -Complex(conductance(1, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 1), -Complex(conductance(1, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 2), -Complex(conductance(1, 2), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 0), -Complex(conductance(2, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 1), -Complex(conductance(2, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 2), -Complex(conductance(2, 2), 0));


		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 0), -Complex(conductance(0, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 1), -Complex(conductance(0, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 2), -Complex(conductance(0, 2), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 0), -Complex(conductance(1, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 1), -Complex(conductance(1, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 2), -Complex(conductance(1, 2), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 0), -Complex(conductance(2, 0), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 1), -Complex(conductance(2, 1), 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 2), -Complex(conductance(2, 2), 0));
	}

	//if (terminalNotGrounded(0))
	//	mSLog->info("Add {} to {}, {}", conductance, matrixNodeIndex(0,0), matrixNodeIndex(0,0));
	//if (terminalNotGrounded(1))
	//	mSLog->info("Add {} to {}, {}", conductance, matrixNodeIndex(1,0), matrixNodeIndex(1,0));
	//if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
	//	mSLog->info("Add {} to {}, {}", -conductance, matrixNodeIndex(0,0), matrixNodeIndex(1,0));
	//	mSLog->info("Add {} to {}, {}", -conductance, matrixNodeIndex(1,0), matrixNodeIndex(0,0));
	//}
}

void DP::Ph3::Resistor::MnaPostStep::execute(Real time, Int timeStepCount) {
	mResistor.mnaUpdateVoltage(**mLeftVector);
	mResistor.mnaUpdateCurrent(**mLeftVector);
}

void DP::Ph3::Resistor::mnaUpdateVoltage(const Matrix& leftVector) {
	// Voltage across component is defined as V1 - V0
	**mIntfVoltage = MatrixComp::Zero(3,1);
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0,0) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1,0));
		(**mIntfVoltage)(1,0) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1,1));
		(**mIntfVoltage)(2,0) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1,2));
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0,0) = (**mIntfVoltage)(0,0) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0,0));
		(**mIntfVoltage)(1,0) = (**mIntfVoltage)(1,0) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0,1));
		(**mIntfVoltage)(2,0) = (**mIntfVoltage)(2,0) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0,2));
	}

	SPDLOG_LOGGER_DEBUG(mSLog, "Voltage A: {} < {}", std::abs((**mIntfVoltage)(0,0)), std::arg((**mIntfVoltage)(0,0)));
}

void DP::Ph3::Resistor::mnaUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = (**mResistance).inverse() * **mIntfVoltage;

	SPDLOG_LOGGER_DEBUG(mSLog, "Current A: {} < {}", std::abs((**mIntfCurrent)(0,0)), std::arg((**mIntfCurrent)(0,0)));
}
