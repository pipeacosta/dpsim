/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/EMT/EMT_Ph3_PiLine.h>

using namespace CPS;

EMT::Ph3::PiLine::PiLine(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	mPhaseType = PhaseType::ABC;
	setVirtualNodeNumber(1);
	setTerminalNumber(2);

	mSLog->info("Create {} {}", this->type(), name);
	mIntfVoltage = Matrix::Zero(3, 1);
	mIntfCurrent = Matrix::Zero(3, 1);

	addAttribute<Matrix>("R_series", &mSeriesRes, Flags::read | Flags::write);
	addAttribute<Matrix>("L_series", &mSeriesInd, Flags::read | Flags::write);
	addAttribute<Matrix>("C_parallel", &mParallelCap, Flags::read | Flags::write);
	addAttribute<Matrix>("G_parallel", &mParallelCond, Flags::read | Flags::write);
	mSLog->flush();
}

SimPowerComp<Real>::Ptr EMT::Ph3::PiLine::clone(String name) {
	auto copy = PiLine::make(name, mLogLevel);
	copy->setParameters(mSeriesRes, mSeriesInd, mParallelCap, mParallelCond);
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
	mParallelCond = (mParallelCond(0, 0) > 0) ? mParallelCond : defaultParallelCond;

	// Static calculation
	Real omega = 2. * PI * frequency;
	MatrixComp impedance = MatrixComp::Zero(3, 3);
	impedance <<
		Complex(mSeriesRes(0, 0), omega * mSeriesInd(0, 0)), Complex(mSeriesRes(0, 1), omega * mSeriesInd(0, 1)), Complex(mSeriesRes(0, 2), omega * mSeriesInd(0, 2)),
		Complex(mSeriesRes(1, 0), omega * mSeriesInd(1, 0)), Complex(mSeriesRes(1, 1), omega * mSeriesInd(1, 1)), Complex(mSeriesRes(1, 2), omega * mSeriesInd(1, 2)),
		Complex(mSeriesRes(2, 0), omega * mSeriesInd(2, 0)), Complex(mSeriesRes(2, 1), omega * mSeriesInd(2, 1)), Complex(mSeriesRes(2, 2), omega * mSeriesInd(2, 2));

	MatrixComp vInitABC = MatrixComp::Zero(3, 1);
	vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(1) - RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
	vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
	vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
	MatrixComp iInit = impedance.inverse() * vInitABC;
	mIntfCurrent = iInit.real();
	mIntfVoltage = vInitABC.real();

	// Initialization of virtual node
	// Initial voltage of phase B,C is set after A
	MatrixComp vInitTerm0 = MatrixComp::Zero(3, 1);
	vInitTerm0(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
	vInitTerm0(1, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_B;
	vInitTerm0(2, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_C;

	mVirtualNodes[0]->setInitialVoltage(PEAK1PH_TO_RMS3PH*(vInitTerm0 + mSeriesRes * iInit));

	// Create series sub components
	mSubSeriesResistor = std::make_shared<EMT::Ph3::Resistor>(mName + "_res", mLogLevel);
	mSubSeriesResistor->setParameters(mSeriesRes);
	mSubSeriesResistor->connect({ mTerminals[0]->node(), mVirtualNodes[0] });
	mSubSeriesResistor->initialize(mFrequencies);
	mSubSeriesResistor->initializeFromNodesAndTerminals(frequency);

	mSubSeriesInductor = std::make_shared<EMT::Ph3::Inductor>(mName + "_ind", mLogLevel);
	mSubSeriesInductor->setParameters(mSeriesInd);
	mSubSeriesInductor->connect({ mVirtualNodes[0], mTerminals[1]->node() });
	mSubSeriesInductor->initialize(mFrequencies);
	mSubSeriesInductor->initializeFromNodesAndTerminals(frequency);

	// Create parallel sub components
	mSubParallelResistor0 = std::make_shared<EMT::Ph3::Resistor>(mName + "_con0", mLogLevel);
	mSubParallelResistor0->setParameters(2. * mParallelCond.inverse());
	mSubParallelResistor0->connect(SimNode::List{ SimNode::GND, mTerminals[0]->node() });
	mSubParallelResistor0->initialize(mFrequencies);
	mSubParallelResistor0->initializeFromNodesAndTerminals(frequency);

	mSubParallelResistor1 = std::make_shared<EMT::Ph3::Resistor>(mName + "_con1", mLogLevel);
	mSubParallelResistor1->setParameters(2. * mParallelCond.inverse());
	mSubParallelResistor1->connect(SimNode::List{ SimNode::GND, mTerminals[1]->node() });
	mSubParallelResistor1->initialize(mFrequencies);
	mSubParallelResistor1->initializeFromNodesAndTerminals(frequency);

	if (mParallelCap(0,0) > 0) {
		mSubParallelCapacitor0 = std::make_shared<EMT::Ph3::Capacitor>(mName + "_cap0", mLogLevel);
		mSubParallelCapacitor0->setParameters(mParallelCap / 2.);
		mSubParallelCapacitor0->connect(SimNode::List{ SimNode::GND, mTerminals[0]->node() });
		mSubParallelCapacitor0->initialize(mFrequencies);
		mSubParallelCapacitor0->initializeFromNodesAndTerminals(frequency);

		mSubParallelCapacitor1 = std::make_shared<EMT::Ph3::Capacitor>(mName + "_cap1", mLogLevel);
		mSubParallelCapacitor1->setParameters(mParallelCap / 2.);
		mSubParallelCapacitor1->connect(SimNode::List{ SimNode::GND, mTerminals[1]->node() });
		mSubParallelCapacitor1->initialize(mFrequencies);
		mSubParallelCapacitor1->initializeFromNodesAndTerminals(frequency);
	}

	mSLog->debug(
		"\n--debug--"
		"\n seriesRes: {:s}"
		"\n seriesInd: {:s}"
		"\n Impedance: {:s}"
		"\n vInit: {:s}"
		"\n iInit: {:s}",
		Logger::matrixToString(mSeriesRes),
		Logger::matrixToString(mSeriesInd),
		Logger::matrixCompToString(impedance),
		Logger::matrixCompToString(vInitABC),
		Logger::matrixCompToString(iInit));

	mSLog->info(
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nTerminal 1 voltage: {:s}"
		"\nVirtual Node 1 voltage: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::matrixToString(mIntfVoltage),
		Logger::matrixToString(mIntfCurrent),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(1)),
		Logger::phasorToString(mVirtualNodes[0]->initialSingleVoltage()));
	mSLog->flush();
}

void EMT::Ph3::PiLine::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();
	MNAInterface::List subComps({ mSubSeriesResistor, mSubSeriesInductor });

	mSubSeriesResistor->mnaInitialize(omega, timeStep, leftVector);
	mSubSeriesInductor->mnaInitialize(omega, timeStep, leftVector);
	mRightVectorStamps.push_back(&mSubSeriesInductor->attribute<Matrix>("right_vector")->get());

	mSubParallelResistor0->mnaInitialize(omega, timeStep, leftVector);
	mSubParallelResistor1->mnaInitialize(omega, timeStep, leftVector);
	subComps.push_back(mSubParallelResistor0);
	subComps.push_back(mSubParallelResistor1);
	
	if (mParallelCap(0,0) > 0) {
		mSubParallelCapacitor0->mnaInitialize(omega, timeStep, leftVector);
		mSubParallelCapacitor1->mnaInitialize(omega, timeStep, leftVector);
		mRightVectorStamps.push_back(&mSubParallelCapacitor0->attribute<Matrix>("right_vector")->get());
		mRightVectorStamps.push_back(&mSubParallelCapacitor1->attribute<Matrix>("right_vector")->get());
		subComps.push_back(mSubParallelCapacitor0);
		subComps.push_back(mSubParallelCapacitor1);
	}
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void EMT::Ph3::PiLine::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	mSubSeriesResistor->mnaApplySystemMatrixStamp(systemMatrix);
	mSubSeriesInductor->mnaApplySystemMatrixStamp(systemMatrix);

	mSubParallelResistor0->mnaApplySystemMatrixStamp(systemMatrix);
	mSubParallelResistor1->mnaApplySystemMatrixStamp(systemMatrix);

	if (mParallelCap(0,0) > 0) {
		mSubParallelCapacitor0->mnaApplySystemMatrixStamp(systemMatrix);
		mSubParallelCapacitor1->mnaApplySystemMatrixStamp(systemMatrix);
	}
}

void EMT::Ph3::PiLine::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	rightVector.setZero();
	for (auto stamp : mRightVectorStamps)
		rightVector += *stamp;
}

void EMT::Ph3::PiLine::mnaAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes){
	// add pre-step dependencies of subcomponents
	mSubSeriesResistor->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	mSubSeriesInductor->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	mSubParallelResistor0->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	mSubParallelResistor1->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	
	if (mParallelCap(0, 0) > 0) {
		mSubParallelCapacitor0->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
		mSubParallelCapacitor1->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	}
	// add pre-step dependencies of component itself
	prevStepDependencies.push_back(attribute("i_intf"));
	prevStepDependencies.push_back(attribute("v_intf"));
	modifiedAttributes.push_back(attribute("right_vector"));
}

void EMT::Ph3::PiLine::mnaPreStep(Real time, Int timeStepCount) {
	// pre-step of subcomponents
	mSubSeriesInductor->mnaPreStep(time, timeStepCount);
	if (mParallelCap(0, 0) > 0) {
		mSubParallelCapacitor0->mnaPreStep(time, timeStepCount);
		mSubParallelCapacitor1->mnaPreStep(time, timeStepCount);
	}
	// pre-step of component itself
	mnaApplyRightSideVectorStamp(mRightVector);
}

void EMT::Ph3::PiLine::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	// add post-step dependencies of subcomponents
	mSubSeriesResistor->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	mSubSeriesInductor->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	if (mParallelCap(0, 0) > 0) {
		mSubParallelCapacitor0->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
		mSubParallelCapacitor1->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	}
	mSubParallelResistor0->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	mSubParallelResistor1->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	// add post-step dependencies of component itself
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(attribute("v_intf"));
	modifiedAttributes.push_back(attribute("i_intf"));
}

void EMT::Ph3::PiLine::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	// post-step of subcomponents
	mSubSeriesResistor->mnaPostStep(time, timeStepCount, leftVector);
	mSubSeriesInductor->mnaPostStep(time, timeStepCount, leftVector);
	if (mParallelCap(0, 0) > 0) {
		mSubParallelCapacitor0->mnaPostStep(time, timeStepCount, leftVector);
		mSubParallelCapacitor1->mnaPostStep(time, timeStepCount, leftVector);
	}
	mSubParallelResistor0->mnaPostStep(time, timeStepCount, leftVector);
	mSubParallelResistor1->mnaPostStep(time, timeStepCount, leftVector);
	// post-step of component itself
	mnaUpdateVoltage(*leftVector);
	mnaUpdateCurrent(*leftVector);
}

void EMT::Ph3::PiLine::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	mIntfVoltage = Matrix::Zero(3, 1);
	if (terminalNotGrounded(1)) {
		mIntfVoltage(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
		mIntfVoltage(1, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 1));
		mIntfVoltage(2, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 2));
	}
	if (terminalNotGrounded(0)) {
		mIntfVoltage(0, 0) = mIntfVoltage(0, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
		mIntfVoltage(1, 0) = mIntfVoltage(1, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
		mIntfVoltage(2, 0) = mIntfVoltage(2, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
	}
}

void EMT::Ph3::PiLine::mnaUpdateCurrent(const Matrix& leftVector) {
	mIntfCurrent = mSubSeriesInductor->intfCurrent();
}

// #### DAE functions ####

void EMT::Ph3::PiLine::daeInitialize(double time, double state[],
	double dstate_dt[], int& offset) {
	// offset: number of component in state, dstate_dt
	// state[offset] = current through voltage source flowing into node matrixNodeIndex(1)
	// dstate_dt[offset] = derivative of current through voltage source  (not used yed) --> set to zero

	updateMatrixNodeIndices();

	// add seriesInductors to state vector
	Matrix inductorCurrent = mSubSeriesInductor->intfCurrent();
	Matrix inductorVoltage = mSubSeriesInductor->intfVoltage();
	state[offset] = inductorCurrent(0,0);
	dstate_dt[offset] = inductorVoltage(0,0)/mSeriesInd(0, 0);
	state[offset+1] = inductorCurrent(1,0);
	dstate_dt[offset+1] = inductorVoltage(1,0)/mSeriesInd(1, 1);
	state[offset+2] = inductorCurrent(2,0);
	dstate_dt[offset+2] = inductorVoltage(2,0)/mSeriesInd(2, 2);


	// Modify initial value of derivative of node voltages
	Matrix capacitor0Current = mSubParallelCapacitor0->intfCurrent();
	Matrix capacitor1Current = mSubParallelCapacitor1->intfCurrent();
	if (mParallelCond(0,0) > 0) {
		if (terminalNotGrounded(1)) {
			dstate_dt[matrixNodeIndex(1, 0)] = capacitor0Current(0,0)/mParallelCap(0,0);
			dstate_dt[matrixNodeIndex(1, 1)] = capacitor0Current(1,0)/mParallelCap(1,1);
			dstate_dt[matrixNodeIndex(1, 2)] = capacitor0Current(2,0)/mParallelCap(2,2);
		}
		if (terminalNotGrounded(0)) {
			dstate_dt[matrixNodeIndex(0, 0)] = capacitor1Current(0,0)/mParallelCap(0,0);
			dstate_dt[matrixNodeIndex(0, 1)] = capacitor1Current(1,0)/mParallelCap(1,1);
			dstate_dt[matrixNodeIndex(0, 2)] = capacitor1Current(2,0)/mParallelCap(2,2);
		}
	}

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nAdded current-phase1 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded current-phase2 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded current-phase3 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded derivative of current-phase1 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current-phase2 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current-phase3 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nSet initial value of derivative of node voltage of capacitor0-phaseA of PiLine '{:s}' ={:f}"
		"\nSet initial value of derivative of node voltage of capacitor0-phaseB of PiLine '{:s}' ={:f}"
		"\nSet initial value of derivative of node voltage of capacitor0-phaseC of PiLine '{:s}' ={:f}"
		"\nSet initial value of derivative of node voltage of capacitor1-phaseA of PiLine '{:s}' ={:f}"
		"\nSet initial value of derivative of node voltage of capacitor1-phaseB of PiLine '{:s}' ={:f}"
		"\nSet initial value of derivative of node voltage of capacitor1-phaseC of PiLine '{:s}' ={:f}"
		"\n--- daeInitialize finished ---",
		this->name(), state[offset],
		this->name(), state[offset+1],
		this->name(), state[offset+2],
		this->name(), dstate_dt[offset],
		this->name(), dstate_dt[offset+1],
		this->name(), dstate_dt[offset+2],
		this->name(), dstate_dt[matrixNodeIndex(0, 0)],
		this->name(), dstate_dt[matrixNodeIndex(0, 1)],
		this->name(), dstate_dt[matrixNodeIndex(0, 2)],
		this->name(), dstate_dt[matrixNodeIndex(1, 0)],
		this->name(), dstate_dt[matrixNodeIndex(1, 1)],
		this->name(), dstate_dt[matrixNodeIndex(1, 2)]
	);
	std::cout << "test\n";
	mSLog->flush();
	offset+=3;
}

void EMT::Ph3::PiLine::daeResidual(double sim_time,
	const double state[], const double dstate_dt[],
	double resid[], std::vector<int>& off) {
	/*
	R*i_rx + L*der(i_rx) = v1 - v2;
    i1 = i_rx + G/2*v1 + C/2*der(v1);
    -i2 = -i_rx + G/2*v2 + C/2*der(v2);
	*/

	int c_offset = off[0]+off[1]; //current offset for component

	resid[c_offset]   = -mSeriesInd(0, 0)*dstate_dt[c_offset]   - state[c_offset]*mSeriesRes(0,0);
	resid[c_offset+1] = -mSeriesInd(1, 1)*dstate_dt[c_offset+1] - state[c_offset]*mSeriesRes(1,1);
	resid[c_offset+2] = -mSeriesInd(2, 2)*dstate_dt[c_offset+2] - state[c_offset]*mSeriesRes(2,2);
	if (terminalNotGrounded(0)) {
		resid[c_offset]   -= state[matrixNodeIndex(0, 0)];
		resid[c_offset+1] -= state[matrixNodeIndex(0, 1)];
		resid[c_offset+2] -= state[matrixNodeIndex(0, 2)];
		resid[matrixNodeIndex(0, 0)] -= state[c_offset] ;
		resid[matrixNodeIndex(0, 1)] -= state[c_offset+1];
		resid[matrixNodeIndex(0, 2)] -= state[c_offset+2];
		resid[matrixNodeIndex(0, 0)] += -state[c_offset] + mParallelCap(0,0)*dstate_dt[matrixNodeIndex(0, 0)]/2 + state[matrixNodeIndex(0, 0)]*mParallelCond(0,0)/2;
		resid[matrixNodeIndex(0, 1)] += -state[c_offset+1] + mParallelCap(1,1)*dstate_dt[matrixNodeIndex(0, 1)]/2 + state[matrixNodeIndex(0, 1)]*mParallelCond(1,1)/2;
		resid[matrixNodeIndex(0, 2)] += -state[c_offset+2] + mParallelCap(2,2)*dstate_dt[matrixNodeIndex(0, 2)]/2 + state[matrixNodeIndex(0, 2)]*mParallelCond(2,2)/2;
	}
	if (terminalNotGrounded(1)) {
		resid[c_offset]   += state[matrixNodeIndex(1, 0)];
		resid[c_offset+1] += state[matrixNodeIndex(1, 1)];
		resid[c_offset+1] += state[matrixNodeIndex(1, 2)];
		resid[matrixNodeIndex(1, 0)] += state[c_offset] + mParallelCap(0,0)*dstate_dt[matrixNodeIndex(1, 0)]/2 + state[matrixNodeIndex(1, 0)]*mParallelCond(0,0)/2;
		resid[matrixNodeIndex(1, 1)] += state[c_offset+1] + mParallelCap(1,1)*dstate_dt[matrixNodeIndex(1, 1)]/2 + state[matrixNodeIndex(1, 1)]*mParallelCond(1,1)/2;
		resid[matrixNodeIndex(1, 2)] += state[c_offset+2] + mParallelCap(2,2)*dstate_dt[matrixNodeIndex(1, 2)]/2 + state[matrixNodeIndex(1, 2)]*mParallelCond(2,2)/2;
	}

	off[1] += 3;

}

void EMT::Ph3::PiLine::daePostStep(const double state[], const double dstate_dt[], int& offset) {
	mIntfVoltage(0,0) = state[matrixNodeIndex(1, 0)] - state[matrixNodeIndex(0, 0)];
	mIntfVoltage(1,0) = state[matrixNodeIndex(1, 1)] - state[matrixNodeIndex(0, 1)];
	mIntfVoltage(2,0) = state[matrixNodeIndex(1, 2)] - state[matrixNodeIndex(0, 2)];
	mIntfCurrent(0,0) = state[offset++];
	mIntfCurrent(1,0) = state[offset++];
	mIntfCurrent(2,0) = state[offset++];
}

