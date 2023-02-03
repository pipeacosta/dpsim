/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_SynchronGenerator3OrderVBR.h>

using namespace CPS;

EMT::Ph3::SynchronGenerator3OrderVBR::SynchronGenerator3OrderVBR
    (const String & uid, const String & name, Logger::Level logLevel)
	: ReducedOrderSynchronGeneratorVBR(uid, name, logLevel),
	mEdq0_t(Attribute<Matrix>::create("Edq0_t", mAttributes)) {

	//
	mSGOrder = SGOrder::SG3Order;

	// model specific variables
	**mEdq0_t = Matrix::Zero(3,1);
	mEhs_vbr = Matrix::Zero(3,1);
}

EMT::Ph3::SynchronGenerator3OrderVBR::SynchronGenerator3OrderVBR
	(const String & name, Logger::Level logLevel)
	: SynchronGenerator3OrderVBR(name, name, logLevel) {
}

void EMT::Ph3::SynchronGenerator3OrderVBR::specificInitialization() {
	// initial voltage behind the transient reactance in the dq0 reference frame
	(**mEdq0_t)(1,0) = (**mVdq0)(1,0) + (**mIdq0)(0,0) * mLd_t;

	mSLog->info(
		"\n--- Model specific initialization  ---"
		"\nInitial Eq_t (per unit): {:f}"
		"\n--- Model specific initialization finished ---",
		(**mEdq0_t)(1,0)
	);
	mSLog->flush();
}

void EMT::Ph3::SynchronGenerator3OrderVBR::stepInPerUnit() {

	if (mSimTime>0.0) {
		// calculate Eq_t at t=k
		(**mEdq0_t)(1,0) = (**mIdq0)(0,0) * mLd_t + (**mVdq0)(1,0);
	}

	// get transformation matrix
	mAbcToDq0 = get_parkTransformMatrix();
	mDq0ToAbc = get_inverseParkTransformMatrix();

	// calculate resistance matrix at t=k+1
	calculateResistanceMatrix();

	// VBR history voltage
	mEhs_vbr(0,0) = 0.0;
	mEhs_vbr(1,0) = mAq_t * (**mIdq0)(0,0) + mBq_t * (**mEdq0_t)(1,0) + mDq_t * mEf_prev + mDq_t * (**mEf);
	mEhs_vbr(2,0) = 0.0;

	// convert Edq_t into the abc reference frame
	mEvbr = mDq0ToAbc * mEhs_vbr * mBase_V;
}

