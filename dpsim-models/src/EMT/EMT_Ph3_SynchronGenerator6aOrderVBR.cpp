/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_SynchronGenerator6aOrderVBR.h>

using namespace CPS;

EMT::Ph3::SynchronGenerator6aOrderVBR::SynchronGenerator6aOrderVBR
    (const String & uid, const String & name, Logger::Level logLevel)
	: ReducedOrderSynchronGeneratorVBR(uid, name, logLevel),
	mEdq0_t(Attribute<Matrix>::create("Edq0_t", mAttributes)),
	mEdq0_s(Attribute<Matrix>::create("Edq0_s", mAttributes))  {

	//
	mSGOrder = SGOrder::SG6aOrder;
	
	// model specific variables
	**mEdq0_t = Matrix::Zero(3,1);
	**mEdq0_s = Matrix::Zero(3,1);
	mEh_t = Matrix::Zero(3,1);
	mEh_s = Matrix::Zero(3,1);
}

EMT::Ph3::SynchronGenerator6aOrderVBR::SynchronGenerator6aOrderVBR
	(const String & name, Logger::Level logLevel)
	: SynchronGenerator6aOrderVBR(name, name, logLevel) {
}

void EMT::Ph3::SynchronGenerator6aOrderVBR::specificInitialization() {
	// initial voltage behind the transient reactance in the dq reference frame
	(**mEdq0_t)(0,0) = (mLq - mLq_t - mYq) * (**mIdq0)(1,0);
	(**mEdq0_t)(1,0) = (1 - mTaa / mTd0_t) * **mEf - (mLd - mLd_t - mYd) * (**mIdq0)(0,0);

	// initial dq behind the subtransient reactance in the dq reference frame
	(**mEdq0_s)(0,0) = (**mVdq0)(0,0) - mLq_s * (**mIdq0)(1,0);
	(**mEdq0_s)(1,0) = (**mVdq0)(1,0) + mLd_s * (**mIdq0)(0,0);

	mSLog->info(
		"\n--- Model specific initialization  ---"
		"\nInitial Ed_t (per unit): {:f}"
		"\nInitial Eq_t (per unit): {:f}"
		"\nInitial Ed_s (per unit): {:f}"
		"\nInitial Eq_s (per unit): {:f}"
		"\n--- Model specific initialization finished ---",

		(**mEdq0_t)(0,0),
		(**mEdq0_t)(1,0),
		(**mEdq0_s)(0,0),
		(**mEdq0_s)(1,0)
	);
	mSLog->flush();
}

void EMT::Ph3::SynchronGenerator6aOrderVBR::stepInPerUnit() {

	if (mSimTime>0.0) {
		// calculate Edq_t at t=k
		(**mEdq0_t)(0,0) = mAd_t * (**mIdq0)(1,0) + mEh_t(0,0);
		(**mEdq0_t)(1,0) = mAq_t * (**mIdq0)(0,0) + mEh_t(1,0);
		(**mEdq0_t)(2,0) = 0.0;

		// calculate Edq_s at t=k
		(**mEdq0_s)(0,0) = -(**mIdq0)(1,0) * mLq_s + (**mVdq0)(0,0);
		(**mEdq0_s)(1,0) = (**mIdq0)(0,0) * mLd_s + (**mVdq0)(1,0);
	}

	// get transformation matrix
	mAbcToDq0 = get_parkTransformMatrix();
	mDq0ToAbc = get_inverseParkTransformMatrix();

	// calculate resistance matrix at t=k+1
	calculateResistanceMatrix();

	// calculate history term behind the transient reactance
	mEh_t(0,0) = mAd_t * (**mIdq0)(1,0) + mBd_t * (**mEdq0_t)(0,0);
	mEh_t(1,0) = mAq_t * (**mIdq0)(0,0) + mBq_t * (**mEdq0_t)(1,0) + mDq_t * (**mEf) + mDq_t * mEf_prev;
	mEh_t(2,0) = 0.0;

	// calculate history term behind the subtransient reactance
	mEh_s(0,0) = mAd_s * (**mIdq0)(1,0) + mBd_s * (**mEdq0_t)(0,0) + mCd_s * (**mEdq0_s)(0,0);
	mEh_s(1,0) = mAq_s * (**mIdq0)(0,0) + mBq_s * (**mEdq0_t)(1,0) + mCq_s * (**mEdq0_s)(1,0) + mDq_s * (**mEf) + mDq_s * mEf_prev;
	mEh_s(2,0) = 0.0;

	// convert Edq_s into the abc reference frame
	mEvbr = mDq0ToAbc * mEh_s * mBase_V;
}

