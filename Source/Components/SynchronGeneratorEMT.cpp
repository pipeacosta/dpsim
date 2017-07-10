#include "SynchronGeneratorEMT.h"

using namespace DPsim;

SynchronGeneratorEMT::SynchronGeneratorEMT(std::string name, int node1, int node2, int node3,
	Real nomPower, Real nomVolt, Real nomFreq, int poleNumber, Real nomFieldCur,
	Real Rs, Real Ll, Real Lmd, Real Lmd0, Real Lmq, Real Lmq0,
	Real Rfd, Real Llfd, Real Rkd, Real Llkd,
	Real Rkq1, Real Llkq1, Real Rkq2, Real Llkq2,
	Real inertia) {

	this->mNode1 = node1 - 1;
	this->mNode2 = node2 - 1;
	this->mNode3 = node3 - 1;

	mNomPower = nomPower;
	mNomVolt = nomVolt;
	mNomFreq = nomFreq;
	mPoleNumber = poleNumber;
	mNomFieldCur = nomFieldCur;

	// base stator values
	mBase_V_RMS = mNomVolt / sqrt(3);
	mBase_v = mBase_V_RMS * sqrt(2);
	mBase_I_RMS = mNomPower / (3 * mBase_V_RMS);
	mBase_i = mBase_I_RMS * sqrt(2);
	mBase_Z = mBase_v / mBase_i;
	mBase_OmElec = 2 * DPS_PI * mNomFreq;
	mBase_OmMech = mBase_OmElec / (mPoleNumber / 2);
	mBase_L = mBase_Z / mBase_OmElec;
	mBase_Psi = mBase_L * mBase_i;
	mBase_T = mNomPower / mBase_OmMech;

	// steady state per unit initial value
	initWithPerUnitParam(Rs, Ll, Lmd, Lmd0, Lmq, Lmq0, Rfd, Llfd, Rkd, Llkd, Rkq1, Llkq1, Rkq2, Llkq2, inertia);
	
}

void SynchronGeneratorEMT::initWithPerUnitParam(
	Real Rs, Real Ll, Real Lmd, Real Lmd0, Real Lmq, Real Lmq0,
	Real Rfd, Real Llfd, Real Rkd, Real Llkd,
	Real Rkq1, Real Llkq1, Real Rkq2, Real Llkq2,
	Real H) {

	// base rotor values
	mBase_ifd = Lmd * mNomFieldCur;
	mBase_vfd = mNomPower / mBase_ifd;
	mBase_Zfd = mBase_vfd / mBase_ifd;
	mBase_Lfd = mBase_Zfd / mBase_OmElec;
		
	mRs = Rs;
	mLl = Ll;
	mLmd = Lmd;
	mLmd0 = Lmd0;
	mLmq = Lmq;
	mLmq0 = Lmq0;
	mRfd = Rfd;
	mLlfd = Llfd;
	mRkd = Rkd;
	mLlkd = Llkd;
	mRkq1 = Rkq1;
	mLlkq1 = Llkq1;
	mRkq2 = Rkq2;
	mLlkq2 = Llkq2;
	mH = H;
	// Additional inductances according to Krause
	mLaq = 1 / (1 / mLmq + 1 / mLl + 1 / mLlkq1 + 1 / mLlkq2);
	mLad = 1 / (1 / mLmd + 1 / mLl + 1 / mLlkd + 1 / mLlfd);

	// Determinant of Ld (inductance matrix of d axis)
	detLd = (mLmd + mLl)*(-mLlfd*mLlkd - mLlfd*mLmd - mLmd*mLlkd) + mLmd*mLmd*(mLlfd + mLlkd);
	// Determinant of Lq (inductance matrix of q axis)
	detLq = -mLmq*mLlkq2*(mLlkq1 + mLl) - mLl*mLlkq1*(mLlkq2 + mLmq);

}

void SynchronGeneratorEMT::init(Real om, Real dt,
	Real initActivePower, Real initReactivePower, Real initTerminalVolt, Real initVoltAngle) {
	
	// steady state per unit initial value
	initStatesInPerUnit(initActivePower, initReactivePower, initTerminalVolt, initVoltAngle);

	mVa = inverseParkTransform2(mThetaMech, mVd* mBase_v, mVq* mBase_v, mV0* mBase_v)(0);
	mVb = inverseParkTransform2(mThetaMech, mVd* mBase_v, mVq* mBase_v, mV0* mBase_v)(1);
	mVc = inverseParkTransform2(mThetaMech, mVd* mBase_v, mVq* mBase_v, mV0* mBase_v)(2);

	mIa = inverseParkTransform2(mThetaMech, mId* mBase_i, mIq* mBase_i, mI0* mBase_i)(0);
	mIb = inverseParkTransform2(mThetaMech, mId* mBase_i, mIq* mBase_i, mI0* mBase_i)(1);
	mIc = inverseParkTransform2(mThetaMech, mId* mBase_i, mIq* mBase_i, mI0* mBase_i)(2);
}

void SynchronGeneratorEMT::initStatesInPerUnit(Real initActivePower, Real initReactivePower,
	Real initTerminalVolt, Real initVoltAngle) {

	double init_P = initActivePower / mNomPower;
	double init_Q = initReactivePower / mNomPower;
	double init_S = sqrt(pow(init_P, 2.) + pow(init_Q, 2.));
	double init_vt = initTerminalVolt / mBase_v;
	double init_it = init_S / init_vt;

	// power factor
	double init_pf = acos(init_P / init_S);

	// load angle
	double init_delta = atan(((mLmq + mLl) * init_it * cos(init_pf) - mRs * init_it * sin(init_pf)) /
		(init_vt + mRs * init_it * cos(init_pf) + (mLmq + mLl) * init_it * sin(init_pf)));
	double init_delta_deg = init_delta / DPS_PI * 180;

	// dq stator voltages and currents
	double init_vd = init_vt * sin(init_delta);
	double init_vq = init_vt * cos(init_delta);
	double init_id = init_it * sin(init_delta + init_pf);
	double init_iq = init_it * cos(init_delta + init_pf);

	// rotor voltage and current
	double init_ifd = (init_vq + mRs * init_iq + (mLmd + mLl) * init_id) / mLmd;
	double init_vfd = mRfd * init_ifd;

	// flux linkages
	double init_psid = init_vq + mRs * init_iq;
	double init_psiq = -init_vd - mRs * init_id;
	double init_psifd = (mLmd + mLlfd) * init_ifd - mLmd * init_id;
	double init_psid1 = mLmd * (init_ifd - init_id);
	double init_psiq1 = -mLmq * init_iq;
	double init_psiq2 = -mLmq * init_iq;

	// rotor mechanical variables
	double init_Te = init_P + mRs * pow(init_it, 2.);
	mOmMech = 1;

	mVd = init_vd;
	mVq = init_vq;
	mV0 = 0;
	mVfd = init_vfd;
	mVkd = 0;
	mVkq1 = 0;
	mVkq2 = 0;

	mIq = init_iq;
	mId = init_id;
	mI0 = 0;
	mIfd = init_ifd;
	mIkd = 0;
	mIkq1 = 0;
	mIkq2 = 0;

	mPsiq = init_psiq;
	mPsid = init_psid;
	mPsi0 = 0;
	mPsifd = init_psifd;
	mPsikd = init_psid1;
	mPsikq1 = init_psiq1;
	mPsikq2 = init_psiq2;

	// Initialize mechanical angle
	//mThetaMech = initVoltAngle + init_delta;
	mThetaMech = initVoltAngle + init_delta - M_PI/2;
}

void SynchronGeneratorEMT::step(SystemModel& system, Real fieldVoltage, Real mechPower, Real time) {

	stepInPerUnit(system.getOmega(), system.getTimeStep(), fieldVoltage, mechPower, time, system.getNumMethod());
	
	// Update current source accordingly
	if (mNode1 >= 0) {
		system.addRealToRightSideVector(mNode1, mIa);
	}
	if (mNode2 >= 0) {
		system.addRealToRightSideVector(mNode2, mIb);
	}
	if (mNode3 >= 0) {
		system.addRealToRightSideVector(mNode3, mIc);
	}

}

void SynchronGeneratorEMT::stepInPerUnit(Real om, Real dt, Real fieldVoltage, Real mechPower, Real time, NumericalMethod numMethod) {
	
	mVa = (1 / mBase_v) * mVa;
	mVb = (1 / mBase_v) * mVb;
	mVc = (1 / mBase_v) * mVc;

	mIa = (1 / mBase_i) * mIa;
	mIb = (1 / mBase_i) * mIb;
	mIc = (1 / mBase_i) * mIc;

	mVfd = fieldVoltage / mBase_v;
	// TODO calculate effect of changed field voltage

	// dq-transform of interface voltage
	mVd = parkTransform2(mThetaMech, mVa, mVb, mVc)(0);
	mVq = parkTransform2(mThetaMech, mVa, mVb, mVc)(1);
	mV0 = parkTransform2(mThetaMech, mVa, mVb, mVc)(2);

	if (numMethod == NumericalMethod::Euler) {

		mMechPower = mechPower / mNomPower;
		mMechTorque = mMechPower / mOmMech;
	
		mElecTorque = (mPsid*mIq - mPsiq*mId);

		// Euler step forward	
		mOmMech = mOmMech + dt * (1 / (2 * mH) * (mMechTorque - mElecTorque));

	}
	else {

		//Two steps Adams-Bashforth
		if (time < dt) {
			// calculate mechanical states
			mMechPower = mechPower / mNomPower;
			mMechTorque = mMechPower / mOmMech;
			mElecTorque = (mPsid*mIq - mPsiq*mId);

			mOmMech_past = mOmMech;
			mOmMech = mOmMech + dt * (1 / (2 * mH) * (mMechTorque - mElecTorque));
		}
		else {
			// calculate mechanical states
			mMechPower = mechPower / mNomPower;
			mMechTorque = mMechPower / mOmMech;
			mMechTorque_past = mMechPower / mOmMech_past;

			mElecTorque = (mPsid*mIq - mPsiq*mId);
			mElecTorque_past = (mPsid_past*mIq - mPsiq_past*mId);
			mOmMech_past = mOmMech;
			mOmMech = mOmMech + (3. / 2.)*dt* (1 / (2 * mH) * (mMechTorque - mElecTorque)) - (1. / 2.)*dt* (1 / (2 * mH) * (mMechTorque_past - mElecTorque_past));
		}
	}
	

	double dtPsid = mVd + mRs*mId + mPsiq*mOmMech;	
	double dtPsiq = mVq + mRs*mIq - mPsid*mOmMech;	
	double dtPsi0 = mV0 + mRs*mI0;
	double dtPsifd = mVfd - mRfd*mIfd;
	double dtPsikd = -mRkd*mIkd;
	double dtPsikq1 = -mRkq1*mIkq1;
	double dtPsikq2 = -mRkq2*mIkq2;

	//if (dtPsid < 0.000001)
	//	dtPsid = 0;
	//if (dtPsiq < 0.000001)
	//	dtPsiq = 0;
	//if (dtPsi0 < 0.000001)
	//	dtPsi0 = 0;
	//if (dtPsifd < 0.000001)
	//	dtPsifd = 0;
	//if (dtPsikd < 0.000001)
	//	dtPsikd = 0;
	//if (dtPsikq1 < 0.000001)
	//	dtPsikq1 = 0;
	//if (dtPsikq2 < 0.000001)
	//	dtPsikq2 = 0;



	mPsid_past = mPsid;
	mPsiq_past = mPsiq;

	mPsid = mPsid + dt*mBase_OmElec*dtPsid;
	mPsiq = mPsiq + dt*mBase_OmElec*dtPsiq;
	mPsi0 = mPsi0 + dt*mBase_OmElec*dtPsi0;
	mPsifd = mPsifd + dt*mBase_OmElec*dtPsifd;
	mPsikd = mPsikd + dt*mBase_OmElec*dtPsikd;
	mPsikq1 = mPsikq1 + dt*mBase_OmElec*dtPsikq1;
	mPsikq2 = mPsikq2 + dt*mBase_OmElec*dtPsikq2;

	

	//Calculation of currents based on inverse of inductance matrix
	mId_past = mId;
	mIq_past = mIq;
	mId = ((mLlfd*mLlkd + mLmd*(mLlfd + mLlkd))*mPsid - mLmd*mLlkd*mPsifd - mLlfd*mLmd*mPsikd) / detLd;
	mIfd = (mLlkd*mLmd*mPsid - (mLl*mLlkd + mLmd*(mLl + mLlkd))*mPsifd + mLmd*mLl*mPsikd) / detLd;
	mIkd = (mLmd*mLlfd*mPsid + mLmd*mLl*mPsifd - (mLmd*(mLlfd + mLl) + mLl*mLlfd)*mPsikd) / detLd;
	mIq = ((mLlkq1*mLlkq2 + mLmq*(mLlkq1 + mLlkq2))*mPsiq - mLmq*mLlkq2*mPsikq1 - mLmq*mLlkq1*mPsikq2) / detLq;
	mIkq1 = (mLmq*mLlkq2*mPsiq - (mLmq*(mLlkq2 + mLl) + mLl*mLlkq2)*mPsikq1 + mLmq*mLl*mPsikq2) / detLq;
	mIkq2 = (mLmq*mLlkq1*mPsiq + mLmq*mLl*mPsikq1 - (mLmq*(mLlkq1 + mLl) + mLl*mLlkq1)*mPsikq2) / detLq;
	mI0 = -mPsi0 / mLl;

	// Update mechanical rotor angle with respect to electrical angle
	mThetaMech = mThetaMech + dt * (mOmMech * mBase_OmMech);


	mIa = mBase_i * inverseParkTransform2(mThetaMech, mId, mIq, mI0)(0);
	mIb = mBase_i * inverseParkTransform2(mThetaMech, mId, mIq, mI0)(1);
	mIc = mBase_i * inverseParkTransform2(mThetaMech, mId, mIq, mI0)(2);

	mCurrents2 << mIq,
		mId,
		mI0,
		mIkq1,
		mIkq2,
		mIfd,
		mIkd;

	mVoltages2 << mVq,
		mVd,
		mV0,
		mVkq1,
		mVkq2,
		mVfd,
		mVkd;

	mFluxes2 << mVq,
		mPsid,
		mPsi0,
		mPsikq1,
		mPsikq2,
		mPsifd,
		mPsikd;
}

void SynchronGeneratorEMT::postStep(SystemModel& system) {
	if (mNode1 >= 0) {
		mVa = system.getRealFromLeftSideVector(mNode1);
	}
	else {
		mVa = 0;
	}
	if (mNode2 >= 0) {
		mVb = system.getRealFromLeftSideVector(mNode2);
	}
	else {
		mVb = 0;
	}
	if (mNode3 >= 0) {
		mVc = system.getRealFromLeftSideVector(mNode3);
	}
	else {
		mVc = 0;
	}
}


DPSMatrix SynchronGeneratorEMT::parkTransform2(Real theta, double a, double b, double c) {
	
	DPSMatrix dq0vector(3, 1);

	// Park transform according to Kundur
	double d, q, zero;

	d = 2. / 3. * cos(theta) * a + 2. / 3. * cos(theta - 2. * M_PI / 3.)*b + 2. / 3. * cos(theta + 2. * M_PI / 3.)*c;
	q = -2. / 3. * sin(theta)*a - 2. / 3. * sin(theta - 2. * M_PI / 3.)*b - 2. / 3. * sin(theta + 2. * M_PI / 3.)*c;
	zero = 1. / 3. * a, 1. / 3. * b, 1. / 3. * c;

	dq0vector << d,
		q,
		0;

	return dq0vector;
}


DPSMatrix SynchronGeneratorEMT::inverseParkTransform2(Real theta, double d, double q, double zero) {
	
	DPSMatrix abcVector(3, 1);

	double a, b, c;

		// Park transform according to Kundur
	a = cos(theta)*d - sin(theta)*q + 1.*zero;
	b = cos(theta - 2. * M_PI / 3.)*d - sin(theta - 2. * M_PI / 3.)*q + 1.*zero;
	c=	cos(theta + 2. * M_PI / 3.)*d - sin(theta + 2. * M_PI / 3.)*q + 1.*zero;

	abcVector << a,
		b,
		c;

	return abcVector;
}