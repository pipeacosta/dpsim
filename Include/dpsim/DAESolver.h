/*********************************************************************************
* @file
* @author Markus Mirz <mmirz@eonerc.rwth-aachen.de>
* @copyright 2017-2018, Institute for Automation of Complex Power Systems, EONERC
*
* CPowerSystems
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

#pragma once 
#include <iostream>
#include <vector>
#include <list>
#include <dpsim/Solver.h>
#include <cps/SystemTopology.h>
//#include <dpsim/Logger.h>
#include <ida/ida.h>
#include <ida/ida_direct.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>

//#define NVECTOR_DATA(vec) NV_DATA_S (vec) // Returns pointer to the first element of array vec

using namespace DPsim ;
using namespace CPS;
	/// Solver class which uses Differential Algebraic Equation(DAE) systems
	class DAESolver : public Solver{
	protected:
		// General simulation parameters
		/// Local copy of the SystemTopology 
        static SystemTopology DAESys;
		/// Offsets vector for adding new equations to the residual vector
        static std::vector<int> offsets;
		/// Constant time step
		Real mTimestep;
		/// Number of equations in problem
		int NEQ;
        ///
        PowerComponent<Real>::List mComponents;
        ///
        static Node<Real>::List mNodes;

		//IDA simulation variables
		/// Memory block allocated by IDA
		void *mem = NULL;
		/// Vector of problem variables
        N_Vector state = NULL;
		///  Derivates of the state vector with respect to time
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
	

	public:
		/// Create solve object with given parameters
        DAESolver(String name,  SystemTopology system, Real dt, Real t0);
        /// Deallocate all memory
        ~DAESolver();
		/// Initialize Components & Nodes with inital values
        void initialize(Real t0);
		/// Residual Function of entire System
        static int DAE_residualFunction(realtype ttime, N_Vector state, N_Vector dstate_dt, N_Vector resid, void *user_data);
		/// Solve system for the current time
		Real step(Real time);
	};

