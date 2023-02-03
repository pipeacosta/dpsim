/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/


#include <dpsim-models/DP/DP_Ph3_Capacitor.h>


using namespace CPS;
using namespace CPS::DP::Ph3;

DP::Ph3::Capacitor::Capacitor(String uid, String name, Logger::Level logLevel)
	: Base::Ph3::Capacitor(mAttributes), SimPowerComp<Complex>(uid, name, logLevel) {
	mPhaseType = PhaseType::ABC;
	setTerminalNumber(2);
	mEquivCurrent = MatrixComp::Zero(3,1);
	**mIntfVoltage = MatrixComp::Zero(3,1);
	**mIntfCurrent = MatrixComp::Zero(3,1);
}

SimPowerComp<Complex>::Ptr DP::Ph3::Capacitor::clone(String name) {
	auto copy = Capacitor::make(name, mLogLevel);
	copy->setParameters(**mCapacitance);
	return copy;
}

void DP::Ph3::Capacitor::initializeFromNodesAndTerminals(Real frequency) {

	Real omega = 2 * PI * frequency;
	MatrixComp susceptance = Matrix::Zero(3, 3);

	susceptance <<
		Complex(0, omega * (**mCapacitance)(0, 0)), Complex(0, omega * (**mCapacitance)(0, 1)), Complex(0, omega * (**mCapacitance)(0, 2)),
		Complex(0, omega * (**mCapacitance)(1, 0)), Complex(0, omega * (**mCapacitance)(1, 1)), Complex(0, omega * (**mCapacitance)(1, 2)),
		Complex(0, omega * (**mCapacitance)(2, 0)), Complex(0, omega * (**mCapacitance)(2, 1)), Complex(0, omega * (**mCapacitance)(2, 2));


	// IntfVoltage initialization for each phase
	(**mIntfVoltage)(0, 0) = initialSingleVoltage(1) - initialSingleVoltage(0);
	Real voltMag = Math::abs((**mIntfVoltage)(0, 0));
	Real voltPhase = Math::phase((**mIntfVoltage)(0, 0));
	(**mIntfVoltage)(1, 0) = Complex(
		voltMag*cos(voltPhase - 2. / 3.*M_PI),
		voltMag*sin(voltPhase - 2. / 3.*M_PI));
	(**mIntfVoltage)(2, 0) = Complex(
		voltMag*cos(voltPhase + 2. / 3.*M_PI),
		voltMag*sin(voltPhase + 2. / 3.*M_PI));

	**mIntfCurrent = susceptance * **mIntfVoltage;

	mSLog->info( "\n--- Initialize from power flow ---" );
				// << "Impedance: " << impedance << std::endl
				// << "Voltage across: " << std::abs((**mIntfVoltage)(0,0))
				// << "<" << Math::phaseDeg((**mIntfVoltage)(0,0)) << std::endl
				// << "Current: " << std::abs((**mIntfCurrent)(0,0))
				// << "<" << Math::phaseDeg((**mIntfCurrent)(0,0)) << std::endl
				// << "Terminal 0 voltage: " << std::abs(initialSingleVoltage(0))
				// << "<" << Math::phaseDeg(initialSingleVoltage(0)) << std::endl
				// << "Terminal 1 voltage: " << std::abs(initialSingleVoltage(1))
				// << "<" << Math::phaseDeg(initialSingleVoltage(1)) << std::endl
				// << "--- Power flow initialization finished ---" << std::endl;
}


void DP::Ph3::Capacitor::initVars(Real omega, Real timeStep) {
	Matrix a = timeStep / 2 * (**mCapacitance).inverse();
	Real b = timeStep * omega / 2.;

	Matrix equivCondReal = a.inverse();
	Matrix equivCondImag = b * equivCondReal;
	mEquivCond = Matrix::Zero(3, 3);
	mEquivCond <<
		Complex(equivCondReal(0, 0), equivCondImag(0, 0)), Complex(equivCondReal(0, 1), equivCondImag(0, 1)), Complex(equivCondReal(0, 2), equivCondImag(0, 2)),
		Complex(equivCondReal(1, 0), equivCondImag(1, 0)), Complex(equivCondReal(1, 1), equivCondImag(1, 1)), Complex(equivCondReal(1, 2), equivCondImag(1, 2)),
		Complex(equivCondReal(2, 0), equivCondImag(2, 0)), Complex(equivCondReal(2, 1), equivCondImag(2, 1)), Complex(equivCondReal(2, 2), equivCondImag(2, 2));

	// since equivCondReal == a.inverse()
	// and
	/*Matrix mPrevVoltCoeffReal = a.inverse();
	Matrix mPrevVoltCoeffImag = -b * a.inverse();*/
	Matrix mPrevVoltCoeffReal = equivCondReal;
	Matrix mPrevVoltCoeffImag = -b * equivCondReal;

	mPrevVoltCoeff = Matrix::Zero(3, 3);
	mPrevVoltCoeff <<
		Complex(mPrevVoltCoeffReal(0, 0), mPrevVoltCoeffImag(0, 0)), Complex(mPrevVoltCoeffReal(0, 1), mPrevVoltCoeffImag(0, 1)), Complex(mPrevVoltCoeffReal(0, 2), mPrevVoltCoeffImag(0, 2)),
		Complex(mPrevVoltCoeffReal(1, 0), mPrevVoltCoeffImag(1, 0)), Complex(mPrevVoltCoeffReal(1, 1), mPrevVoltCoeffImag(1, 1)), Complex(mPrevVoltCoeffReal(1, 2), mPrevVoltCoeffImag(1, 2)),
		Complex(mPrevVoltCoeffReal(2, 0), mPrevVoltCoeffImag(2, 0)), Complex(mPrevVoltCoeffReal(2, 1), mPrevVoltCoeffImag(2, 1)), Complex(mPrevVoltCoeffReal(2, 2), mPrevVoltCoeffImag(2, 2));


	mEquivCurrent = -mPrevVoltCoeff * **mIntfVoltage - **mIntfCurrent;
}

void DP::Ph3::Capacitor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	updateMatrixNodeIndices();
	initVars(omega, timeStep);
	//Matrix equivCondReal = 2.0 * mCapacitance / timeStep;
	//Matrix equivCondImag = omega * mCapacitance;
	//mEquivCond <<
	//	Complex(equivCondReal(0, 0), equivCondImag(0, 0)),
	//	Complex(equivCondReal(1, 0), equivCondImag(1, 0)),
	//	Complex(equivCondReal(2, 0), equivCondImag(2, 0));

	// TODO: something is wrong here -- from Ph1_Capacitor
	/*Matrix prevVoltCoeffReal = 2.0 * mCapacitance / timeStep;
	Matrix prevVoltCoeffImag = - omega * mCapacitance;
	mPrevVoltCoeff = Matrix::Zero(3, 1);
	mPrevVoltCoeff <<
		Complex(prevVoltCoeffReal(0, 0), prevVoltCoeffImag(0, 0)),
		Complex(prevVoltCoeffReal(1, 0), prevVoltCoeffImag(1, 0)),
		Complex(prevVoltCoeffReal(2, 0), prevVoltCoeffImag(2, 0));

	mEquivCurrent = -**mIntfCurrent + -mPrevVoltCoeff.cwiseProduct( **mIntfVoltage);*/
	// no need to update current now
	//**mIntfCurrent = mEquivCond.cwiseProduct(**mIntfVoltage) + mEquivCurrent;

	// mLog.info() << "\n--- MNA Initialization ---" << std::endl
	// 			<< "Initial voltage " << Math::abs((**mIntfVoltage)(0,0))
	// 			<< "<" << Math::phaseDeg((**mIntfVoltage)(0,0)) << std::endl
	// 			<< "Initial current " << Math::abs((**mIntfCurrent)(0,0))
	// 			<< "<" << Math::phaseDeg((**mIntfCurrent)(0,0)) << std::endl
	// 			<< "--- MNA initialization finished ---" << std::endl;

	**mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

void DP::Ph3::Capacitor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {

	if (terminalNotGrounded(0)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0), mEquivCond(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 1), mEquivCond(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 2), mEquivCond(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 0), mEquivCond(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 1), mEquivCond(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 2), mEquivCond(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 0), mEquivCond(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 1), mEquivCond(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 2), mEquivCond(2, 2));
	}
	if (terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 0), mEquivCond(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 1), mEquivCond(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 2), mEquivCond(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 0), mEquivCond(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 1), mEquivCond(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 2), mEquivCond(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 0), mEquivCond(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 1), mEquivCond(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 2), mEquivCond(2, 2));
	}
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 0), -mEquivCond(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 1), -mEquivCond(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 2), -mEquivCond(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 0), -mEquivCond(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 1), -mEquivCond(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 2), -mEquivCond(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 0), -mEquivCond(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 1), -mEquivCond(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 2), -mEquivCond(2, 2));

		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 0), -mEquivCond(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 1), -mEquivCond(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 2), -mEquivCond(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 0), -mEquivCond(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 1), -mEquivCond(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 2), -mEquivCond(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 0), -mEquivCond(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 1), -mEquivCond(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 2), -mEquivCond(2, 2));
	}/*
	mLog.debug() << "\n--- Apply system matrix stamp ---" << std::endl;
	if (terminalNotGrounded(0)) {
		mLog.debug() << "Add " << mEquivCond(0, 0) << " to " << matrixNodeIndex(0, 0) << "," << matrixNodeIndex(0, 0) << std::endl;
		mLog.debug() << "Add " << mEquivCond(1, 0) << " to " << matrixNodeIndex(0, 1) << "," << matrixNodeIndex(0, 1) << std::endl;
		mLog.debug() << "Add " << mEquivCond(2, 0) << " to " << matrixNodeIndex(0, 2) << "," << matrixNodeIndex(0, 2) << std::endl;
	}
	if (terminalNotGrounded(1)) {
		mLog.debug() << "Add " << mEquivCond(0, 0) << " to " << matrixNodeIndex(1, 0) << "," << matrixNodeIndex(1, 0) << std::endl;
		mLog.debug() << "Add " << mEquivCond(0, 1) << " to " << matrixNodeIndex(1, 1) << "," << matrixNodeIndex(1, 1) << std::endl;
		mLog.debug() << "Add " << mEquivCond(0, 2) << " to " << matrixNodeIndex(1, 2) << "," << matrixNodeIndex(1, 2) << std::endl;
	}
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		mLog.debug() << "Add " << -mEquivCond(0, 0) << " to " << matrixNodeIndex(0, 0) << "," << matrixNodeIndex(1, 0) << std::endl
			<< "Add " << -mEquivCond(0, 0) << " to " << matrixNodeIndex(1, 0) << "," << matrixNodeIndex(0, 0) << std::endl;
		mLog.debug() << "Add " << -mEquivCond(1, 0) << " to " << matrixNodeIndex(0, 1) << "," << matrixNodeIndex(1, 1) << std::endl
			<< "Add " << -mEquivCond(1, 0) << " to " << matrixNodeIndex(1, 1) << "," << matrixNodeIndex(0, 1) << std::endl;
		mLog.debug() << "Add " << -mEquivCond(2, 0) << " to " << matrixNodeIndex(0, 2) << "," << matrixNodeIndex(1, 2) << std::endl
			<< "Add " << -mEquivCond(2, 0) << " to " << matrixNodeIndex(1, 2) << "," << matrixNodeIndex(0, 2) << std::endl;
	}*/
}

void DP::Ph3::Capacitor::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	//mCureqr = mCurrr + mGcr * mDeltavr + mGci * mDeltavi;
	//mCureqi = mCurri + mGcr * mDeltavi - mGci * mDeltavr;

	mEquivCurrent = -**mIntfCurrent + -mPrevVoltCoeff * **mIntfVoltage;

	if (terminalNotGrounded(0)) {
		Math::setVectorElement(rightVector, matrixNodeIndex(0, 0), mEquivCurrent(0, 0));
		Math::setVectorElement(rightVector, matrixNodeIndex(0, 1), mEquivCurrent(1, 0));
		Math::setVectorElement(rightVector, matrixNodeIndex(0, 2), mEquivCurrent(2, 0));
	}
	if (terminalNotGrounded(1)) {
		Math::setVectorElement(rightVector, matrixNodeIndex(1, 0), -mEquivCurrent(0, 0));
		Math::setVectorElement(rightVector, matrixNodeIndex(1, 1), -mEquivCurrent(1, 0));
		Math::setVectorElement(rightVector, matrixNodeIndex(1, 2), -mEquivCurrent(2, 0));
	}
}

void DP::Ph3::Capacitor::MnaPreStep::execute(Real time, Int timeStepCount) {
	mCapacitor.mnaApplyRightSideVectorStamp(**mCapacitor.mRightVector);
}

void DP::Ph3::Capacitor::MnaPostStep::execute(Real time, Int timeStepCount) {
	mCapacitor.mnaUpdateVoltage(**mLeftVector);
	mCapacitor.mnaUpdateCurrent(**mLeftVector);
}

void DP::Ph3::Capacitor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	**mIntfVoltage = Matrix::Zero(3, 1);
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0, 0) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1, 0));
		(**mIntfVoltage)(1, 0) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1, 1));
		(**mIntfVoltage)(2, 0) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1, 2));
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0, 0) = (**mIntfVoltage)(0, 0) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0, 0));
		(**mIntfVoltage)(1, 0) = (**mIntfVoltage)(1, 0) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0, 1));
		(**mIntfVoltage)(2, 0) = (**mIntfVoltage)(2, 0) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0, 2));
	}
}

void DP::Ph3::Capacitor::mnaUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = mEquivCond * **mIntfVoltage + mEquivCurrent;
}
