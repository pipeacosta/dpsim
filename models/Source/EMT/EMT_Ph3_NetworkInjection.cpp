/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <cps/EMT/EMT_Ph3_NetworkInjection.h>

using namespace CPS;

EMT::Ph3::NetworkInjection::NetworkInjection(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	mPhaseType = PhaseType::ABC;
	setVirtualNodeNumber(0);
	setTerminalNumber(1);
	mIntfVoltage = Matrix::Zero(3, 1);
	mIntfCurrent = Matrix::Zero(3, 1);

	mSLog->info("Create {} {}", this->type(), name);

	// Create electrical sub components
	mSubVoltageSource = std::make_shared<EMT::Ph3::VoltageSource>(mName + "_vs", mLogLevel);
	mSubComponents.push_back(mSubVoltageSource);
	mSLog->info("Electrical subcomponents: ");
	for (auto subcomp: mSubComponents)
		mSLog->info("- {}", subcomp->name());

	addAttributeRef<MatrixComp>("V_ref", mSubVoltageSource->attribute<MatrixComp>("V_ref"), Flags::read | Flags::write);
	addAttributeRef<Real>("f_src", mSubVoltageSource->attribute<Real>("f_src"), Flags::read | Flags::write);
}

SimPowerComp<Real>::Ptr EMT::Ph3::NetworkInjection::clone(String name) {
	auto copy = NetworkInjection::make(name, mLogLevel);
	copy->setParameters(attribute<MatrixComp>("V_ref")->get());
	return copy;
}

void EMT::Ph3::NetworkInjection::setParameters(MatrixComp voltageRef, Real srcFreq) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(voltageRef, srcFreq);

	mSLog->info("\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]", 
				Logger::matrixCompToString(voltageRef),
				Logger::realToString(srcFreq));
}

void EMT::Ph3::NetworkInjection::initializeFromNodesAndTerminals(Real frequency) {
	// Connect electrical subcomponents
	mSubVoltageSource->connect({ SimNode::GND, node(0) });
	
	// Initialize electrical subcomponents
	for (auto subcomp: mSubComponents) {
		subcomp->initialize(mFrequencies);
		subcomp->initializeFromNodesAndTerminals(frequency);
	}
}

// #### MNA functions ####

void EMT::Ph3::NetworkInjection::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	// initialize electrical subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaInitialize(omega, timeStep, leftVector);

	// collect right side vectors of subcomponents
	mRightVectorStamps.push_back(&mSubVoltageSource->attribute<Matrix>("right_vector")->get());

	// collect tasks
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));

	mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void EMT::Ph3::NetworkInjection::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaApplySystemMatrixStamp(systemMatrix);
}

void EMT::Ph3::NetworkInjection::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	rightVector.setZero();
	for (auto stamp : mRightVectorStamps)
		rightVector += *stamp;

	mSLog->debug("Right Side Vector: {:s}",
				Logger::matrixToString(rightVector));
}


void EMT::Ph3::NetworkInjection::mnaAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
	// add pre-step dependencies of subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	// add pre-step dependencies of component itself
	prevStepDependencies.push_back(attribute("i_intf"));
	prevStepDependencies.push_back(attribute("v_intf"));
	modifiedAttributes.push_back(attribute("right_vector"));
}

void EMT::Ph3::NetworkInjection::mnaPreStep(Real time, Int timeStepCount) {
	// pre-step of subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaPreStep(time, timeStepCount);
	// pre-step of component itself
	mnaApplyRightSideVectorStamp(mRightVector);
}

void EMT::Ph3::NetworkInjection::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	// add post-step dependencies of subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	// add post-step dependencies of component itself
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(attribute("v_intf"));
	modifiedAttributes.push_back(attribute("i_intf"));
}

void EMT::Ph3::NetworkInjection::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	// post-step of subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaPostStep(time, timeStepCount, leftVector);
	// post-step of component itself
	mnaUpdateCurrent(*leftVector);
	mnaUpdateVoltage(*leftVector);
}

void EMT::Ph3::NetworkInjection::mnaUpdateVoltage(const Matrix& leftVector) {
	mIntfVoltage = mSubVoltageSource->attribute<Matrix>("v_intf")->get();
}

void EMT::Ph3::NetworkInjection::mnaUpdateCurrent(const Matrix& leftVector) {
	mIntfCurrent = mSubVoltageSource->attribute<Matrix>("i_intf")->get();
}


// #### DAE functions ####

void EMT::Ph3::NetworkInjection::daeInitialize(double time, double state[],
	double dstate_dt[], int& offset) {
	//current is positive when flows out of the network injection

	updateMatrixNodeIndices();
	this->updateVoltage(time);
	state[offset] = mIntfCurrent(0,0);
	dstate_dt[offset] = 0.0;
	state[offset+1] = mIntfCurrent(1,0);
	dstate_dt[offset+1] = 0.0;
	state[offset+2] = mIntfCurrent(2,0);
	dstate_dt[offset+2] = 0.0;

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nInitial time={:f}s"
		"\nAdded current phase1 of NetworkInjection '{:s}' to state vector, initial value={:f}A"
		"\nAdded current phase2 of NetworkInjection '{:s}' to state vector, initial value={:f}A"
		"\nAdded current phase3 of NetworkInjection '{:s}' to state vector, initial value={:f}A"
		"\nAdded derivative of current phase1 of NetworkInjection '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current phase2 of NetworkInjection '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current phase3 of NetworkInjection '{:s}' to derivative state vector, initial value={:f}"
		"\n--- daeInitialize finished ---",
		time,
		this->name(), state[offset],
		this->name(), state[offset+1],
		this->name(), state[offset+2],
		this->name(), dstate_dt[offset],
		this->name(), dstate_dt[offset+1],
		this->name(), dstate_dt[offset+2]
	);
	mSLog->flush();
	offset+=3;
}

void EMT::Ph3::NetworkInjection::daeResidual(double sim_time,
	const double state[], const double dstate_dt[],
	double resid[], std::vector<int>& off) {

	this->updateVoltage(sim_time);
	int c_offset = off[0]+off[1]; //current offset for component
	
	int pos_node1 = matrixNodeIndex(0, 0);
	int pos_node2 = matrixNodeIndex(0, 1);
	int pos_node3 = matrixNodeIndex(0, 2);
	resid[c_offset]   = state[pos_node1] - mIntfVoltage(0,0);
	resid[c_offset+1] = state[pos_node2] - mIntfVoltage(1,0);
	resid[c_offset+2] = state[pos_node3] - mIntfVoltage(2,0);

	// update nodal equations
	resid[pos_node1] -= state[c_offset];
	resid[pos_node2] -= state[c_offset+1];
	resid[pos_node3] -= state[c_offset+2];

	mSLog->debug(
		"\n\n--- NetworkInjection name: {:s} - SimTime= {:f} ---"
		"\nresid[c_offset]   = state[pos_node1] - mIntfVoltage(0,0) = {:f} - {:f} = {:f}"
		"\nresid[c_offset+1] = state[pos_node2] - mIntfVoltage(1,0) = {:f} - {:f} = {:f}"
		"\nresid[c_offset+2] = state[pos_node3] - mIntfVoltage(2,0) = {:f} - {:f} = {:f}"

		"\nupdate nodal equations:"
		"\nresid[pos_node1] -= state[c_offset]   --> resid[pos_node1] -= {:f}"
		"\nresid[pos_node2] -= state[c_offset+1] --> resid[pos_node2] -= {:f}"
		"\nresid[pos_node3] -= state[c_offset+2] --> resid[pos_node3] -= {:f}",

		this->name(), sim_time,
		state[pos_node1], mIntfVoltage(0,0), resid[c_offset],
		state[pos_node2], mIntfVoltage(1,0), resid[c_offset+1],
		state[pos_node3], mIntfVoltage(2,0), resid[c_offset+2],
		state[c_offset],
		state[c_offset+1],
		state[c_offset+2] 
	);
	//mSLog->flush();

	off[1]+=3;
}

void EMT::Ph3::NetworkInjection::daePostStep(double Nexttime, const double state[], 
	const double dstate_dt[], int& offset) {
	//this->updateVoltage(Nexttime);
	mIntfCurrent(0,0) = state[offset++];
	mIntfCurrent(1,0) = state[offset++];
	mIntfCurrent(2,0) = state[offset++];
}
