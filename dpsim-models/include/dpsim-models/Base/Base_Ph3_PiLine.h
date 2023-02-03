/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/Definitions.h>
#include <dpsim-models/AttributeList.h>

namespace CPS {
namespace Base {
namespace Ph3 {
	class PiLine {

	public:
		/// Resistance along the line [ohms]
		const Attribute<Matrix>::Ptr mSeriesRes;
		/// Inductance along the line [H]
		const Attribute<Matrix>::Ptr mSeriesInd;
		/// Capacitance in parallel to the line [F]
		const Attribute<Matrix>::Ptr mParallelCap;
		/// Conductance in parallel to the line [S]
		const Attribute<Matrix>::Ptr mParallelCond;

		explicit PiLine(CPS::AttributeBase::Map &attributeList) :
			mSeriesRes(Attribute<Matrix>::create("R_series", attributeList)),
			mSeriesInd(Attribute<Matrix>::create("L_series", attributeList)),
			mParallelCap(Attribute<Matrix>::create("C_parallel", attributeList)),
			mParallelCond(Attribute<Matrix>::create("G_parallel", attributeList)) { };

		///
		void setParameters(Matrix seriesResistance, Matrix seriesInductance,
			Matrix parallelCapacitance = Matrix::Zero(3,3), Matrix parallelConductance = Matrix::Zero(3,3)) const {
			**mSeriesRes = seriesResistance;
			**mSeriesInd = seriesInductance;
			**mParallelCond = parallelConductance;
			**mParallelCap = parallelCapacitance;
		}
	};
}
}
}
