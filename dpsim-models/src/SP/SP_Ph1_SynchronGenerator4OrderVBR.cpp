/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/SP/SP_Ph1_SynchronGenerator4OrderVBR.h>
#include <dpsim-models/Signal/TurbineGovernor.h>

using namespace CPS;

SP::Ph1::SynchronGenerator4OrderVBR::SynchronGenerator4OrderVBR
    (const String & uid, const String & name, Logger::Level logLevel)
	: ReducedOrderSynchronGeneratorVBR(uid, name, logLevel),
	mEdq_t(Attribute<Matrix>::create("Edq_t", mAttributes)) {

	//
	mSGOrder = SGOrder::SG4Order;

	// model specific variables
	**mEdq_t = Matrix::Zero(2,1);
	mEh_vbr = Matrix::Zero(2,1);
}

SP::Ph1::SynchronGenerator4OrderVBR::SynchronGenerator4OrderVBR
	(const String & name, Logger::Level logLevel)
	: SynchronGenerator4OrderVBR(name, name, logLevel) {
}

void SP::Ph1::SynchronGenerator4OrderVBR::specificInitialization() {

	// initial voltage behind the transient reactance in the dq reference frame
	(**mEdq_t)(0,0) = (**mVdq)(0,0) - (**mIdq)(1,0) * mLq_t;
	(**mEdq_t)(1,0) = (**mVdq)(1,0) + (**mIdq)(0,0) * mLd_t;

	mSLog->info(
		"\n--- Model specific initialization  ---"
		"\nInitial Ed_t (per unit): {:f}"
		"\nInitial Eq_t (per unit): {:f}"
		"\n--- Model specific initialization finished ---",

		(**mEdq_t)(0,0),
		(**mEdq_t)(1,0)
	);
	mSLog->flush();
}

void SP::Ph1::SynchronGenerator4OrderVBR::stepInPerUnit() {
	
	if (mSimTime>0.0) {
		// calculate Edq_t at t=k
		(**mEdq_t)(0,0) = -(**mIdq)(1,0) * mLq_t + (**mVdq)(0,0);
		(**mEdq_t)(1,0) = (**mIdq)(0,0) * mLd_t + (**mVdq)(1,0);
	}

	// get transformation matrix
	mDqToComplexA = get_DqToComplexATransformMatrix();
	mComplexAToDq = mDqToComplexA.transpose();

	// calculate resistance matrix at t=k+1
	calculateResistanceMatrix();

	// VBR history voltage
	mEh_vbr(0,0) = mAd_t * (**mIdq)(1,0) + mBd_t * (**mEdq_t)(0,0);
	mEh_vbr(1,0) = mAq_t * (**mIdq)(0,0) + mBq_t * (**mEdq_t)(1,0) + mDq_t * mEf_prev + mDq_t * (**mEf);
	
	// convert Edq_t into the abc reference frame
	mEh_vbr = mDqToComplexA * mEh_vbr;
	mEvbr = Complex(mEh_vbr(0,0), mEh_vbr(1,0)) * mBase_V_RMS;
}