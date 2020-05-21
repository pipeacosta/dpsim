/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim/DAESolver.h>

using namespace DPsim;
using namespace CPS;

static int check_retval(void *returnvalue, const char *funcname, int opt);

template <typename VarType>
DAESolver<VarType>::DAESolver(String name, 
    CPS::SystemTopology system, Real dt, 
    Real t0, CPS::Logger::Level logLevel) :
	Solver(name, logLevel), mSystem(system),
	mTimestep(dt) {

    // Defines offset vector of the residual which is composed as follows:
    // mOffset[0] = # nodal voltage equations
    // mOffset[1] = # of components equations (1 for inductor, cap and voltagesource)
    mOffsets.push_back(0);
    mOffsets.push_back(0);

    // Set initial values of all required variables and create IDA solver environment
    mSLog->info("-- Process system components");
    this->initializeComponents();

    mSLog->info("-- Process system nodes");
    for (auto baseNode : mSystem.mNodes) {
        // Add nodes to the list and ignore ground nodes.
        if (!baseNode->isGround()) {
            auto node = std::dynamic_pointer_cast<CPS::SimNode<VarType>>(baseNode);
            if (!node) {
                throw CPS::Exception(); 
            }
            mNodes.push_back(node);
            mNEQ +=1;
            mSLog->info("Added node {:s}", node->name());;
        }
    }

    UInt matrixNodeIndexIdx = 0;
    for (UInt idx = 0; idx < mNodes.size(); idx++) {
        mNodes[idx]->setMatrixNodeIndex(0, matrixNodeIndexIdx);
        matrixNodeIndexIdx++;
        if (mNodes[idx]->phaseType() == PhaseType::ABC) {
            mNodes[idx]->setMatrixNodeIndex(1, matrixNodeIndexIdx);
            matrixNodeIndexIdx++;
            mNodes[idx]->setMatrixNodeIndex(2, matrixNodeIndexIdx);
            matrixNodeIndexIdx++;
        }
    }

    mSLog->info("");
    mSLog->flush();
    initialize(t0);
}

template <typename VarType>
void DAESolver<VarType>::initializeComponents() {
    mSLog->info("-- Initialize components from power flow");
    for(IdentifiedObject::Ptr comp : mSystem.mComponents) {
        // Initialize componentes and add componente eqs. to state vector
        auto emtComp = std::dynamic_pointer_cast<SimPowerComp<VarType>>(comp);
        if (!emtComp) {
            throw CPS::Exception();
        }
        // Set initial values of all components
        emtComp->initializeFromPowerflow(mSystem.mSystemFrequency);
        
        auto daeComp = std::dynamic_pointer_cast<DAEInterface>(comp);
        if (!daeComp) {
            throw CPS::Exception();
        }
        mDAEComponents.push_back(daeComp);

        // Register residual functions of components
        mResidualFunctions.push_back(
                [daeComp](double sim_time, const double state[], const double dstate_dt[], 
                            double resid[], std::vector<int> &off) {
                                daeComp->daeResidual(sim_time, state, dstate_dt, resid, off);
                            });
        
        mNEQ += daeComp->get_numberOfStateVariables();
        mSLog->info("Added {:s} '{:s}' to simulation.", comp->type(), comp->name());
        mSLog->flush();
    }
}

template <typename VarType>
void DAESolver<VarType>::initialize(Real t0) {
    mSLog->info("---- Start dae initialization ----");
    mSLog->info("Number of Eqn.: {}", mNEQ);
    
    mSimTime = t0;
    int counter = 0;
    realtype *sval = NULL, *s_dtval = NULL;
    
    // creates and allocates memory for the state and dstate N Vectors
    state = N_VNew_Serial(mNEQ);
    if(check_retval((void *)state, "N_VNew_Serial", 0)) 
        throw CPS::Exception();
    dstate_dt = N_VNew_Serial(mNEQ);
    if(check_retval((void *)dstate_dt, "N_VNew_Serial", 0)) 
        throw CPS::Exception();

    // capturing a returned array/pointer
    sval  = N_VGetArrayPointer(state);
    s_dtval = N_VGetArrayPointer_Serial(dstate_dt);
    for (auto node : mNodes) {
        // Initialize nodal voltages of state vector
        Real tempVolt = std::real(node->initialSingleVoltage(node->phaseType()));
        sval[counter] = tempVolt;
        s_dtval[counter] = 0;
        mSLog->info("Init voltage node {:s} = {}V", node->name(), sval[counter]);
        counter++;
    }
    for (auto daeComp : mDAEComponents) {
        // Initialize component voltages of state vector
        daeComp->daeInitialize(t0, sval, s_dtval, counter);
    }

    for (int i = 0; i < (mNEQ); i++) {
        // Set initial values for state derivative for now all equal to 0
        s_dtval[i] = 0; // TODO: add derivative calculation
        mSLog->info("derivative {}: {}", i, s_dtval[i]);
    }

    // Set relative tolerance and absolute error
    rtol = RCONST(1.0e-10); // Set relative tolerance
    abstol = RCONST(1.0e-4); // Set absolute error

    // creates the IDA solver memory block
    mem = IDACreate();
    if (mem == NULL)  {
        std::cout << "Error: IDACreate() returns NULL" << std::endl;
        throw CPS::Exception();
    }
    mSLog->info("Creates the IDA solver memory block Ok");
    
    mSLog->info("Define Userdata");
    // This passes the solver instance as the user_data argument to the residual functions
    int ret = IDASetUserData(mem, this);
	if (check_retval(&ret, "IDASetUserData", 1)) {
		throw CPS::Exception();
	}

    mSLog->info("Call IDAInit");
    ret = IDAInit(mem, &DAESolver::residualFunctionWrapper, t0, state, dstate_dt);
    if (check_retval(&ret, "IDAInit", 1)) {
		throw CPS::Exception();
	}

    mSLog->info("Call IDATolerances");
    ret = IDASStolerances(mem, rtol, abstol);
    if (check_retval(&ret, "IDASStolerances", 1)) {
    	throw CPS::Exception();
	}

    /*
    mSLog->info("Call IDASetStopTime");
    ret = IDASetStopTime(mem, mTimestep);
    if (check_retval(&ret, "IDASetStopTime", 1)) {
    	throw CPS::Exception();
	}
    */

    mSLog->info("Call IDA Solver Stuff");
    // Allocate and connect Matrix A and solver LS to IDA
    A = SUNDenseMatrix(mNEQ, mNEQ);
    LS = SUNDenseLinearSolver(state, A);
    ret = IDADlsSetLinearSolver(mem, LS, A);

    //Optional IDA input functions
    //ret = IDASetMaxNumSteps(mem, -1);  //Max. number of timesteps until tout (-1 = unlimited)
    //ret = IDASetMaxConvFails(mem, 100); //Max. number of convergence failures at one step

    mSLog->info("--- Finished initialization --- \n");
    mSLog->flush();
}

template <typename VarType>
int DAESolver<VarType>::residualFunctionWrapper(realtype time, 
    N_Vector state, N_Vector dstate_dt, N_Vector resid, void *user_data)
{
    DAESolver *self = reinterpret_cast<DAESolver *>(user_data);
    return self->residualFunction(time, state, dstate_dt, resid);
}

template <typename VarType>
int DAESolver<VarType>::residualFunction(realtype time, 
    N_Vector state, N_Vector dstate_dt, N_Vector resid)
{
    mOffsets[0] = mNodes.size(); // Reset Offset of nodes
    mOffsets[1] = 0;             // Reset Offset of componentes

    //reset residual functions of nodes
    realtype *residual = NULL;
    residual  = N_VGetArrayPointer(resid);
    for (int i=0; i<mOffsets[0]; i++) {
        residual[i] = 0;
    }

    // Call all registered component residual functions
    for (auto resFn : mResidualFunctions) {
        resFn(mSimTime, NV_DATA_S(state), NV_DATA_S(dstate_dt), NV_DATA_S(resid), mOffsets);
    }

    // If successful; positive value if recoverable error, negative if fatal error
    // TODO: Error handling
    return 0;
}

template <typename VarType>
Real DAESolver<VarType>::step(Real time) {
    realtype NextTime = (realtype) time+mTimestep;

    int ret = IDASolve(mem, NextTime, &tret, state, dstate_dt, IDA_NORMAL);  // TODO: find alternative to IDA_NORMAL
    if (ret != IDA_SUCCESS) {
        mSLog->info("Ida Error: {}", ret);
        void(IDAGetNumSteps(mem, &interalSteps));
        void(IDAGetNumResEvals(mem, &resEval));
        mSLog->info("Interal steps: {}", interalSteps);
        mSLog->info("Res Eval: {}", resEval);
        throw CPS::Exception();
    } 

    //update node voltages
    realtype *sval = NULL;
    realtype *dstate_val = NULL;
    sval  = N_VGetArrayPointer(state);
    dstate_val  = N_VGetArrayPointer(dstate_dt);
    mOffsets[0] = 0;             // Reset Offset of nodes
    for (auto node : mNodes) {
        node->setVoltage(sval[mOffsets[0]]);
        mOffsets[0] += 1;
    }

    //update components
    mOffsets[1] = mOffsets[0];      // Reset Offset of componentes
    for (auto comp : mDAEComponents) { 
        comp->daePostStep(sval, dstate_val, mOffsets[1]);
    }

    mSimTime = time+mTimestep;
    return mSimTime;
}

template <typename VarType>
Task::List DAESolver<VarType>::getTasks() {
    Task::List l;
    l.push_back(std::make_shared<DAESolver<VarType>::SolveStep>(*this));
    l.push_back(std::make_shared<DAESolver<VarType>::LogTask>(*this));
    return l;
}

template <typename VarType>
DAESolver<VarType>::~DAESolver() {
    // Releasing all memory allocated by IDA
    IDAFree(&mem);
    N_VDestroy(state);
    N_VDestroy(dstate_dt);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
}

static int check_retval(void *returnvalue, const char *funcname, int opt) {
    int *retval;
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        std::cout << "\nSUNDIALS_ERROR: " << funcname << "() failed - returned NULL pointer" << std::endl;
        return(1);
    } else if (opt == 1) {
        /* Check if retval < 0 */
        retval = (int *) returnvalue;
        if (*retval < 0) {
            std::cout << "\nSUNDIALS_ERROR: " << funcname << "() failed with retval = " << *retval << std::endl;
            return(1);
        }
    } else if (opt == 2 && returnvalue == NULL) {
        /* Check if function returned NULL pointer - no memory allocated */
         std::cout << "\nMEMORY_ERROR: " << funcname << "() failed - returned NULL pointer" << std::endl;
        return(1);
    }

    return (0);
}

template class DPsim::DAESolver<Real>;
template class DPsim::DAESolver<Complex>;