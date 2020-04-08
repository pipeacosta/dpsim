/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <iostream>
#include <vector>
#include <list>

#include <dpsim/Solver.h>
#include <cps/Solver/DAEInterface.h>
#include <dpsim/Scheduler.h>
#include <cps/SimPowerComp.h>
#include <cps/Logger.h>
#include <dpsim/DataLogger.h>
#include <cps/AttributeList.h>

#include <ida/ida.h>
#include <ida/ida_direct.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>


namespace DPsim {

	/// Solver class which uses Differential Algebraic Equation(DAE) systems
	class DAESolver : public Solver {
	protected:
		// General simulation parameters
        CPS::SystemTopology mSystem;
		/// Offsets vector for adding new equations to the residual vector
		std::vector<Int> mOffsets;
		/// Constant time step
		Real mTimestep;
		/// Number of equations in problem
		Int mNEQ;
		/// Components of the Problem
        CPS::IdentifiedObject::List mComponents;
		CPS::DAEInterface::List mComponents2;		//change name

		/// Nodes of the Problem
        CPS::SimNode<Real>::List mNodes;
		// CPS::SimNode<Complex>::List mNodes;

		// IDA simulation variables
		/// Memory block allocated by IDA
		void *mem = NULL;
		/// Vector of problem variables
		N_Vector state = NULL;
		/// Derivates of the state vector with respect to time
		N_Vector dstate_dt = NULL;
		/// Time IDA reached while solving
		realtype tret;
		/// Scalar absolute tolerance
		realtype abstol;
		/// Relative tolerance
		realtype rtol;
		/// Template Jacobian Matrix
		SUNMatrix A = NULL;
		/// Linear solver object
		SUNLinearSolver LS = NULL;
        long int interalSteps = 0;
        long int resEval=0;
        std::vector<CPS::DAEInterface::ResFn> mResidualFunctions;

		// #### Attributes related to logging ####
		///residual vector logger
		std::shared_ptr<DataLogger> mResidualLog;

		/// Residual Function of entire System
		static int residualFunctionWrapper(realtype ttime, N_Vector state, N_Vector dstate_dt, N_Vector resid, void *user_data);
		int residualFunction(realtype ttime, N_Vector state, N_Vector dstate_dt, N_Vector resid);

	public:
		/// Create solve object with given parameters
        DAESolver(String name, 
			CPS::SystemTopology system, Real dt, Real t0, 
			CPS::Logger::Level logLevel = CPS::Logger::Level::info);
		/// Deallocate all memory
		~DAESolver();
		/// Initialize Components & Nodes with inital values
		void initialize(Real t0);
		/// Solve system for the current time
		Real step(Real time);

		CPS::Task::List getTasks();

		class SolveStep : public CPS::Task {
		public:
			SolveStep(DAESolver& solver) :
				Task(solver.mName + ".SolveStep"), mSolver(solver) {
				mModifiedAttributes.push_back(Scheduler::external);
				//for (auto node : solver.mNodes) {
				//	mModifiedAttributes.push_back(node->attribute("v"));
				//}
			}
			void execute(Real time, Int timeStepCount) {
    			mSolver.step(time);
			}
		private:
			DAESolver& mSolver;
		};

		///
		class LogTask : public CPS::Task {
		public:
			LogTask(DAESolver& solver) :
				Task(solver.mName + ".Log"), mSolver(solver) {
				// mAttributeDependencies.push_back(solver.attribute("left_vector"));
				mModifiedAttributes.push_back(Scheduler::external);
			}

			void execute(Real time, Int timeStepCount) {
				mSolver.log(time);
			}

		private:
			DAESolver& mSolver;
		};

		/// Solution vector of unknown quantities
		Matrix mResidualVector;
		Matrix& ResidualVector() { return mResidualVector; }	///Getter

		/// Log residual vector values for each simulation step
		void log(Real time);
	};
}
