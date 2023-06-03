/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>

using namespace CPS;

EMT::Ph3::PiLine::PiLine(String uid, String name, Logger::Level logLevel)
	: Base::Ph3::PiLine(mAttributes), CompositePowerComp<Real>(uid, name, true, true, logLevel) {
	mPhaseType = PhaseType::ABC;
	setVirtualNodeNumber(1);
	setTerminalNumber(2);

	SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);
	**mIntfVoltage = Matrix::Zero(3, 1);
	**mIntfCurrent = Matrix::Zero(3, 1);

	mSLog->flush();
}

/// DEPRECATED: Delete method
SimPowerComp<Real>::Ptr EMT::Ph3::PiLine::clone(String name) {
	auto copy = PiLine::make(name, mLogLevel);
	copy->setParameters(**mSeriesRes, **mSeriesInd, **mParallelCap, **mParallelCond);
	return copy;
}

void EMT::Ph3::PiLine::initializeFromNodesAndTerminals(Real frequency) {

	// By default there is always a small conductance to ground to
	// avoid problems with floating nodes.
	Matrix defaultParallelCond = Matrix::Zero(3, 3);
	defaultParallelCond <<
		1e-6, 0, 0,
		0, 1e-6, 0,
		0, 0, 1e-6;
	**mParallelCond = ((**mParallelCond)(0, 0) > 0) ? **mParallelCond : defaultParallelCond;

	// Static calculation
	Real omega = 2. * PI * frequency;
	MatrixComp impedance = MatrixComp::Zero(3, 3);
	impedance <<
		Complex((**mSeriesRes)(0, 0), omega * (**mSeriesInd)(0, 0)), Complex((**mSeriesRes)(0, 1), omega * (**mSeriesInd)(0, 1)), Complex((**mSeriesRes)(0, 2), omega * (**mSeriesInd)(0, 2)),
		Complex((**mSeriesRes)(1, 0), omega * (**mSeriesInd)(1, 0)), Complex((**mSeriesRes)(1, 1), omega * (**mSeriesInd)(1, 1)), Complex((**mSeriesRes)(1, 2), omega * (**mSeriesInd)(1, 2)),
		Complex((**mSeriesRes)(2, 0), omega * (**mSeriesInd)(2, 0)), Complex((**mSeriesRes)(2, 1), omega * (**mSeriesInd)(2, 1)), Complex((**mSeriesRes)(2, 2), omega * (**mSeriesInd)(2, 2));

	MatrixComp vInitABC = MatrixComp::Zero(3, 1);
	vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(1) - RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
	vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
	vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
	MatrixComp iInit = impedance.inverse() * vInitABC;
	**mIntfCurrent = iInit.real();
	**mIntfVoltage = vInitABC.real();

	// Initialization of virtual node
	// Initial voltage of phase B,C is set after A
	MatrixComp vInitTerm0 = MatrixComp::Zero(3, 1);
	vInitTerm0(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
	vInitTerm0(1, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_B;
	vInitTerm0(2, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_C;

	mVirtualNodes[0]->setInitialVoltage(PEAK1PH_TO_RMS3PH*(vInitTerm0 + **mSeriesRes * iInit));

	// Create series sub components
	mSubSeriesResistor = std::make_shared<EMT::Ph3::Resistor>(**mName + "_res", mLogLevel);
	mSubSeriesResistor->setParameters(**mSeriesRes);
	mSubSeriesResistor->connect({ mTerminals[0]->node(), mVirtualNodes[0] });
	mSubSeriesResistor->initialize(mFrequencies);
	mSubSeriesResistor->initializeFromNodesAndTerminals(frequency);
	addMNASubComponent(mSubSeriesResistor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

	mSubSeriesInductor = std::make_shared<EMT::Ph3::Inductor>(**mName + "_ind", mLogLevel);
	mSubSeriesInductor->setParameters(**mSeriesInd);
	mSubSeriesInductor->connect({ mVirtualNodes[0], mTerminals[1]->node() });
	mSubSeriesInductor->initialize(mFrequencies);
	mSubSeriesInductor->initializeFromNodesAndTerminals(frequency);
	addMNASubComponent(mSubSeriesInductor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

	// Create parallel sub components
	mSubParallelResistor0 = std::make_shared<EMT::Ph3::Resistor>(**mName + "_con0", mLogLevel);
	mSubParallelResistor0->setParameters(2. * (**mParallelCond).inverse());
	mSubParallelResistor0->connect(SimNode::List{ SimNode::GND, mTerminals[0]->node() });
	mSubParallelResistor0->initialize(mFrequencies);
	mSubParallelResistor0->initializeFromNodesAndTerminals(frequency);
	addMNASubComponent(mSubParallelResistor0, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

	mSubParallelResistor1 = std::make_shared<EMT::Ph3::Resistor>(**mName + "_con1", mLogLevel);
	mSubParallelResistor1->setParameters(2. * (**mParallelCond).inverse());
	mSubParallelResistor1->connect(SimNode::List{ SimNode::GND, mTerminals[1]->node() });
	mSubParallelResistor1->initialize(mFrequencies);
	mSubParallelResistor1->initializeFromNodesAndTerminals(frequency);
	addMNASubComponent(mSubParallelResistor1, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);

	if ((**mParallelCap)(0,0) > 0) {
		mSubParallelCapacitor0 = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_cap0", mLogLevel);
		mSubParallelCapacitor0->setParameters(**mParallelCap / 2.);
		mSubParallelCapacitor0->connect(SimNode::List{ SimNode::GND, mTerminals[0]->node() });
		mSubParallelCapacitor0->initialize(mFrequencies);
		mSubParallelCapacitor0->initializeFromNodesAndTerminals(frequency);
		addMNASubComponent(mSubParallelCapacitor0, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

		mSubParallelCapacitor1 = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_cap1", mLogLevel);
		mSubParallelCapacitor1->setParameters(**mParallelCap / 2.);
		mSubParallelCapacitor1->connect(SimNode::List{ SimNode::GND, mTerminals[1]->node() });
		mSubParallelCapacitor1->initialize(mFrequencies);
		mSubParallelCapacitor1->initializeFromNodesAndTerminals(frequency);
		addMNASubComponent(mSubParallelCapacitor1, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
	}

	SPDLOG_LOGGER_DEBUG(mSLog, 
		"\n--debug--"
		"\n seriesRes: {:s}"
		"\n seriesInd: {:s}"
		"\n Impedance: {:s}"
		"\n vInit: {:s}"
		"\n iInit: {:s}",
		Logger::matrixToString(**mSeriesRes),
		Logger::matrixToString(**mSeriesInd),
		Logger::matrixCompToString(impedance),
		Logger::matrixCompToString(vInitABC),
		Logger::matrixCompToString(iInit));

	SPDLOG_LOGGER_INFO(mSLog, 
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nTerminal 1 voltage: {:s}"
		"\nVirtual Node 1 voltage: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::matrixToString(**mIntfVoltage),
		Logger::matrixToString(**mIntfCurrent),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(1)),
		Logger::phasorToString(mVirtualNodes[0]->initialSingleVoltage()));
	mSLog->flush();
}

void EMT::Ph3::PiLine::mnaParentAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes){
	prevStepDependencies.push_back(mIntfCurrent);
	prevStepDependencies.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::PiLine::mnaParentPreStep(Real time, Int timeStepCount) {
	mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::PiLine::mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mIntfCurrent);
}

void EMT::Ph3::PiLine::mnaParentPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	mnaCompUpdateVoltage(**leftVector);
	mnaCompUpdateCurrent(**leftVector);
}

void EMT::Ph3::PiLine::mnaCompUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	**mIntfVoltage = Matrix::Zero(3, 1);
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
		(**mIntfVoltage)(1, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 1));
		(**mIntfVoltage)(2, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 2));
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0, 0) = (**mIntfVoltage)(0, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
		(**mIntfVoltage)(1, 0) = (**mIntfVoltage)(1, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
		(**mIntfVoltage)(2, 0) = (**mIntfVoltage)(2, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
	}
}

void EMT::Ph3::PiLine::mnaCompUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = mSubSeriesInductor->intfCurrent();
}
