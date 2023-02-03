/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/SP/SP_Ph1_NetworkInjection.h>

using namespace CPS;


SP::Ph1::NetworkInjection::NetworkInjection(String uid, String name,
    Logger::Level logLevel) : CompositePowerComp<Complex>(uid, name, true, true, logLevel),
	mVoltageRef(Attribute<Complex>::createDynamic("V_ref", mAttributes)),
	mSrcFreq(Attribute<Real>::createDynamic("f_src", mAttributes)),
	mVoltageSetPoint(Attribute<Real>::create("V_set", mAttributes)),
	mVoltageSetPointPerUnit(Attribute<Real>::create("V_set_pu", mAttributes, 1.0)),
	mActivePowerInjection(Attribute<Real>::create("p_inj", mAttributes)),
	mReactivePowerInjection(Attribute<Real>::create("q_inj", mAttributes)) {

	mSLog->info("Create {} of type {}", **mName, this->type());
	mSLog->flush();
	**mIntfVoltage = MatrixComp::Zero(1, 1);
	**mIntfCurrent = MatrixComp::Zero(1, 1);
	setVirtualNodeNumber(0);
	setTerminalNumber(1);

	// Create electrical sub components
	mSubVoltageSource = std::make_shared<SP::Ph1::VoltageSource>(**mName + "_vs", mLogLevel);
	addMNASubComponent(mSubVoltageSource, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
	mSLog->info("Electrical subcomponents: ");
	for (auto subcomp: mSubComponents)
		mSLog->info("- {}", subcomp->name());

	// MNA attributes
	mSubVoltageSource->mVoltageRef->setReference(mVoltageRef);
	mSubVoltageSource->mSrcFreq->setReference(mSrcFreq);
}

// #### Powerflow section ####

void SP::Ph1::NetworkInjection::setParameters(Real voltageSetPoint) {
	**mVoltageSetPoint = voltageSetPoint;

	mSLog->info("Voltage Set-Point ={} [V]", **mVoltageSetPoint);
	mSLog->flush();

	mParametersSet = true;
}

void SP::Ph1::NetworkInjection::setParameters(Complex initialPhasor, Real freqStart, Real rocof, Real timeStart, Real duration, bool smoothRamp) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(initialPhasor, freqStart, rocof, timeStart, duration, smoothRamp);

	mSLog->info("\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::phasorToString(initialPhasor),
				Logger::realToString(freqStart));
}

void SP::Ph1::NetworkInjection::setParameters(Complex initialPhasor, Real modulationFrequency, Real modulationAmplitude, Real baseFrequency /*= 0.0*/, bool zigzag /*= false*/) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(initialPhasor, modulationFrequency, modulationAmplitude, baseFrequency, zigzag);

	mSLog->info("\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::phasorToString(initialPhasor),
				Logger::realToString(baseFrequency));
}

void SP::Ph1::NetworkInjection::setBaseVoltage(Real baseVoltage) {
    mBaseVoltage = baseVoltage;
}

void SP::Ph1::NetworkInjection::calculatePerUnitParameters(Real baseApparentPower, Real baseOmega) {
    mSLog->info("#### Calculate Per Unit Parameters for {}", **mName);
	mSLog->info("Base Voltage={} [V]", mBaseVoltage);

    **mVoltageSetPointPerUnit = **mVoltageSetPoint / mBaseVoltage;

	mSLog->info("Voltage Set-Point ={} [pu]", **mVoltageSetPointPerUnit);
	mSLog->flush();
}

void SP::Ph1::NetworkInjection::modifyPowerFlowBusType(PowerflowBusType powerflowBusType) {
	mPowerflowBusType = powerflowBusType;
}

void SP::Ph1::NetworkInjection::updatePowerInjection(Complex powerInj) {
	**mActivePowerInjection = powerInj.real();
	**mReactivePowerInjection = powerInj.imag();
}

// #### MNA Section ####

void SP::Ph1::NetworkInjection::setParameters(Complex voltageRef, Real srcFreq) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(voltageRef, srcFreq);

	mSLog->info("\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::phasorToString(voltageRef),
				Logger::realToString(srcFreq));
}

SimPowerComp<Complex>::Ptr SP::Ph1::NetworkInjection::clone(String name) {
	auto copy = NetworkInjection::make(name, mLogLevel);
	copy->setParameters(**mVoltageRef);
	return copy;
}

void SP::Ph1::NetworkInjection::initializeFromNodesAndTerminals(Real frequency) {
	// Connect electrical subcomponents
	mSubVoltageSource->connect({ SimNode::GND, node(0) });

	// Initialize electrical subcomponents
	for (auto subcomp: mSubComponents) {
		subcomp->initialize(mFrequencies);
		subcomp->initializeFromNodesAndTerminals(frequency);
	}
}

// #### MNA functions ####
void SP::Ph1::NetworkInjection::mnaParentApplyRightSideVectorStamp(Matrix& rightVector) {
	mSLog->debug("Right Side Vector: {:s}",
				Logger::matrixToString(rightVector));
}

void SP::Ph1::NetworkInjection::mnaParentAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
	prevStepDependencies.push_back(mIntfCurrent);
	prevStepDependencies.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mRightVector);
}

void SP::Ph1::NetworkInjection::mnaParentPreStep(Real time, Int timeStepCount) {
	mnaApplyRightSideVectorStamp(**mRightVector);
}

void SP::Ph1::NetworkInjection::mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mIntfCurrent);
}

void SP::Ph1::NetworkInjection::mnaParentPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	mnaUpdateCurrent(**leftVector);
	mnaUpdateVoltage(**leftVector);
}

void SP::Ph1::NetworkInjection::mnaUpdateVoltage(const Matrix& leftVector) {
	**mIntfVoltage = **mSubVoltageSource->mIntfVoltage;
}

void SP::Ph1::NetworkInjection::mnaUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = **mSubVoltageSource->mIntfCurrent;
}

void SP::Ph1::NetworkInjection::daeResidual(double ttime, const double state[], const double dstate_dt[], double resid[], std::vector<int>& off) {
	/* new state vector definintion:
		state[0]=node0_voltage
		state[1]=node1_voltage
		....
		state[n]=noden_voltage
		state[n+1]=component0_voltage
		state[n+2]=component0_inductance (not yet implemented)
		...
		state[m-1]=componentm_voltage
		state[m]=componentm_inductance
	*/

	int Pos1 = matrixNodeIndex(0);
	int Pos2 = matrixNodeIndex(1);
	int c_offset = off[0] + off[1]; //current offset for component
	int n_offset_1 = c_offset + Pos1 + 1;// current offset for first nodal equation
	int n_offset_2 = c_offset + Pos2 + 1;// current offset for second nodal equation
	resid[c_offset] = (state[Pos2] - state[Pos1]) - state[c_offset]; // Voltage equation for Resistor
	//resid[++c_offset] = ; //TODO : add inductance equation
	resid[n_offset_1] += (**mIntfCurrent)(0, 0).real();
	resid[n_offset_2] += (**mIntfCurrent)(0, 0).real();
	off[1] += 1;
}

Complex SP::Ph1::NetworkInjection::daeInitialize() {
	(**mIntfVoltage)(0,0) = (**mSubVoltageSource->mIntfVoltage)(0,0);
	return (**mSubVoltageSource->mIntfVoltage)(0,0);
}
