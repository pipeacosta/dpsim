/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/Base/Base_ReducedOrderSynchronGenerator.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>

namespace CPS {
namespace DP {
namespace Ph1 {
	/// @brief Base class for DP VBR synchronous generator model single phase
	class ReducedOrderSynchronGeneratorVBR :
		public Base::ReducedOrderSynchronGenerator<Complex>,
		public MNAVariableCompInterface {

	public:
        // Common elements of all VBR models
		/// voltage behind reactance phase a
        Complex mEvbr;
		/// norton equivalent current of mEvbr
		Complex mIvbr;

	protected:
        /// Resistance matrix in dq reference frame
		Matrix mResistanceMatrixDq;
		/// Conductance matrix phase A
		Matrix mConductanceMatrix;
		/// Ka Matrix
		MatrixComp mKa;
		/// Kb Matrix
		MatrixComp mKb;
		/// Kb Matrix
		MatrixComp mKc;
		/// Constant part of Resistance matrix in abc reference frame
		Matrix mResistanceMatrix_const;
		/// phase a equivalent part of mKa
		Complex mKa_1ph;
		/// phase a equivalent part of mKb
		Complex mKb_1ph;
		/// phase a equivalent part of mResistanceMatrix_const
		Complex mR_const_1ph;
		/// Vector to create abc vector from a component
		MatrixComp mShiftVector;
		/// complex conjugate of mShiftVector
		MatrixComp mShiftVectorConj;
		/// Matrix to convert Evbr from dq domain to abc domain (only phase a)
		MatrixComp mKvbr;

        /// Constructor 
        ReducedOrderSynchronGeneratorVBR(const String & uid, const String & name, Logger::Level logLevel);
        ReducedOrderSynchronGeneratorVBR(const String & name, Logger::Level logLevel);
      
	  	// #### General Functions ####
        /// Specific component initialization
        virtual void specificInitialization() override =0;
        ///
		void initializeResistanceMatrix() override;
        ///
        virtual void stepInPerUnit() override =0;
		///
        void calculateConductanceMatrix();
		/// Calculate Ka, Kb and Kvbr
		void calculateAuxiliarVariables();
		///
		Matrix get_parkTransformMatrix() const;        

		// ### MNA Section ###
		void mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector);
		void mnaApplySystemMatrixStamp(Matrix& systemMatrix) override;
        void mnaApplyRightSideVectorStamp(Matrix& rightVector) override;
		void mnaPostStep(const Matrix& leftVector) override;
		
    public:
        virtual ~ReducedOrderSynchronGeneratorVBR()=default;

        /// Mark that parameter changes so that system matrix is updated
		Bool hasParameterChanged() override { return true; };
    };
}
}
}