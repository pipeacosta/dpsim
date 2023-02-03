/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_SeriesSwitch.h>

using namespace CPS;

// !!! TODO: 	Adaptions to use in EMT_Ph3 models phase-to-ground peak variables
// !!! 			with initialization from phase-to-phase RMS variables

EMT::Ph3::SeriesSwitch::SeriesSwitch(String uid, String name, Logger::Level logLevel)
	: Base::Ph1::Switch(mAttributes), SimPowerComp<Real>(uid, name, logLevel) {
	mPhaseType = PhaseType::ABC;
	setTerminalNumber(2);
}

SimPowerComp<Real>::Ptr EMT::Ph3::SeriesSwitch::clone(String name) {
	auto copy = SeriesSwitch::make(name, mLogLevel);
	copy->setParameters(**mOpenResistance, **mClosedResistance);
	return copy;
}

void EMT::Ph3::SeriesSwitch::initializeFromNodesAndTerminals(Real frequency) {

	Real impedance = (**mIsClosed) ? **mClosedResistance : **mOpenResistance;

	Complex phasorA = initialSingleVoltage(1) - initialSingleVoltage(0);
	(**mIntfVoltage)(0,0) = phasorA.real();
	Complex alpha(cos(2./3.*PI), sin(2./3.*PI));
	(**mIntfVoltage)(1, 0) = Complex(phasorA * pow(alpha,2)).real();
	(**mIntfVoltage)(2, 0) = Complex(phasorA * alpha).real();

	**mIntfCurrent = **mIntfVoltage / impedance;

	mSLog->info("\n--- Initialization from powerflow ---"
		"\nVoltage across amplitude and phase: \n{}"
		"\nCurrent amplitude and phase: \n{}"
		"\nTerminal 0 voltage amplitude and phase: \n{}"
		"\nTerminal 1 voltage amplitude and phase: \n{}"
		"\n--- Initialization from powerflow finished ---",
		Logger::phasorMatrixToString(**mIntfVoltage),
		Logger::phasorMatrixToString(**mIntfCurrent),
		Logger::phasorMatrixToString(initialVoltage(0)),
		Logger::phasorMatrixToString(initialVoltage(1)));
}

void EMT::Ph3::SeriesSwitch::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

Bool EMT::Ph3::SeriesSwitch::mnaIsClosed() { return **mIsClosed; }

void EMT::Ph3::SeriesSwitch::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	Real conductance = (**mIsClosed)
		? 1. / **mClosedResistance
		: 1. / **mOpenResistance;

	// Set diagonal entries
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndices(0), matrixNodeIndices(0), conductance);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndices(1), matrixNodeIndices(1), conductance);
	// Set off diagonal entries
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndices(0), matrixNodeIndices(1), -conductance);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndices(1), matrixNodeIndices(0), -conductance);
	}

	if (terminalNotGrounded(0))
		mSLog->info("Add {} to {}, {}", conductance, matrixNodeIndices(0)[0], matrixNodeIndices(0)[0]);
	if (terminalNotGrounded(1))
		mSLog->info("Add {} to {}, {}", conductance, matrixNodeIndices(1)[0], matrixNodeIndices(1)[0]);
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		mSLog->info("Add {} to {}, {}", -conductance, matrixNodeIndices(0)[0], matrixNodeIndices(1)[0]);
		mSLog->info("Add {} to {}, {}", -conductance, matrixNodeIndices(1)[0], matrixNodeIndices(0)[0]);
	}
}

void EMT::Ph3::SeriesSwitch::mnaApplySwitchSystemMatrixStamp(Bool closed, Matrix& systemMatrix, Int freqIdx) {
	Real conductance = (closed)
		? 1. / **mClosedResistance
		: 1. / **mOpenResistance;

	// Set diagonal entries
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndices(0), matrixNodeIndices(0), conductance);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndices(1), matrixNodeIndices(1), conductance);
	// Set off diagonal entries
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndices(0), matrixNodeIndices(1), -conductance);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndices(1), matrixNodeIndices(0), -conductance);
	}

	if (terminalNotGrounded(0))
		mSLog->info("Add {} to {}, {}", conductance, matrixNodeIndices(0)[0], matrixNodeIndices(0)[0]);
	if (terminalNotGrounded(1))
		mSLog->info("Add {} to {}, {}", conductance, matrixNodeIndices(1)[0], matrixNodeIndices(1)[0]);
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		mSLog->info("Add {} to {}, {}", -conductance, matrixNodeIndices(0)[0], matrixNodeIndices(1)[0]);
		mSLog->info("Add {} to {}, {}", -conductance, matrixNodeIndices(1)[0], matrixNodeIndices(0)[0]);
	}
}

void EMT::Ph3::SeriesSwitch::MnaPostStep::execute(Real time, Int timeStepCount) {
	mSwitch.mnaUpdateVoltage(**mLeftVector);
	mSwitch.mnaUpdateCurrent(**mLeftVector);
}

void EMT::Ph3::SeriesSwitch::mnaUpdateVoltage(const Matrix& leftVector) {
	// Voltage across component is defined as V1 - V0
	**mIntfVoltage = Matrix::Zero(3,1);
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1,0));
		(**mIntfVoltage)(1,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1,1));
		(**mIntfVoltage)(2,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1,2));
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0,0) = (**mIntfVoltage)(0,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0,0));
		(**mIntfVoltage)(1,0) = (**mIntfVoltage)(1,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0,1));
		(**mIntfVoltage)(2,0) = (**mIntfVoltage)(2,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0,2));
	}

	SPDLOG_LOGGER_DEBUG(mSLog, "Voltage A: {}", (**mIntfVoltage)(0,0));
}

void EMT::Ph3::SeriesSwitch::mnaUpdateCurrent(const Matrix& leftVector) {
	Real impedance = (**mIsClosed)? **mClosedResistance : **mOpenResistance;
	**mIntfCurrent = **mIntfVoltage / impedance;

	SPDLOG_LOGGER_DEBUG(mSLog, "Current A: {}", (**mIntfCurrent)(0,0));
}
