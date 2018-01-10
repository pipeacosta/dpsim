/** InductorDP (EMT)
 *
 * @author Markus Mirz <mmirz@eonerc.rwth-aachen.de>
 * @copyright 2017, Institute for Automation of Complex Power Systems, EONERC
 *
 * DPsim
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *********************************************************************************/

#include "EMT_Inductor.h"

using namespace DPsim;

Components::EMT::Inductor::Inductor(String name, Int src, Int dest, Real inductance)
	: Base(name, src, dest)
{
	mInductance = inductance;
	attrMap["inductance"] = { Attribute::Real, &mInductance };
}

void Components::EMT::Inductor::applySystemMatrixStamp(SystemModel& system)
{
	mGl = system.getTimeStep() / (2.0 * mInductance);

	if (mNode1 >= 0) {
		system.addRealToSystemMatrix(mNode1, mNode1, mGl);
	}
	if (mNode2 >= 0) {
		system.addRealToSystemMatrix(mNode2, mNode2, mGl);
	}
	if (mNode1 >= 0 && mNode2 >= 0) {
		system.addRealToSystemMatrix(mNode1, mNode2, -mGl);
		system.addRealToSystemMatrix(mNode2, mNode1, -mGl);
	}
}


void Components::EMT::Inductor::initialize(SystemModel& system) {
	mCurr = 0;
	mCureq = 0;
	mDeltav = 0;
}

void Components::EMT::Inductor::step(SystemModel& system, Real time)
{
	// Initialize internal state
	mCureq = mGl * mDeltav + mCurr;

	if (mNode1 >= 0) {
		system.addRealToRightSideVector(mNode1, -mCureq);
	}
	if (mNode2 >= 0) {
		system.addRealToRightSideVector(mNode2, mCureq);
	}
}

void Components::EMT::Inductor::postStep(SystemModel& system)
{
	Real vpos, vneg;

	if (mNode1 >= 0) {
		vpos = system.getRealFromLeftSideVector(mNode1);
	}
	else {
		vpos = 0;
	}
	if (mNode2 >= 0) {
		vneg = system.getRealFromLeftSideVector(mNode2);
	}
	else {
		vneg = 0;
	}
	mDeltav = vpos - vneg;
	mCurr = mGl * mDeltav + mCureq;
}