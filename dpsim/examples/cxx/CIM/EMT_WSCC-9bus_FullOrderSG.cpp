/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <iostream>
#include <list>

#include <DPsim.h>

using namespace DPsim;
using namespace CPS::EMT;

int main(int argc, char *argv[]) {

	// Simulation parameters
	String simName = "EMT_WSCC-9bus_FullOrderSG";
	Real timeStep;
	Real finalTime;

	// Find CIM files
	std::list<fs::path> filenames;
	CommandLineArgs args(argc, argv);
	if (argc <= 1) {
		filenames = Utils::findFiles({
			"WSCC-09_Dyn_Full_DI.xml",
			"WSCC-09_Dyn_Full_EQ.xml",
			"WSCC-09_Dyn_Full_SV.xml",
			"WSCC-09_Dyn_Full_TP.xml"
		}, "/cimdata/WSCC-09_Dyn_Full", "CIMPATH");
		timeStep = 10e-6;
		finalTime = 0.1;
	}
	else {
		filenames = args.positionalPaths();
		timeStep = args.timeStep;
		finalTime = args.duration;
	}

	// ----- POWERFLOW FOR INITIALIZATION -----
	// read original network topology
	String simNamePF = simName + "_PF";
	Logger::setLogDir("logs/" + simNamePF);
    CPS::CIM::Reader reader(simNamePF, Logger::Level::debug, Logger::Level::debug);
    SystemTopology systemPF = reader.loadCIM(60, filenames, Domain::SP, PhaseType::Single, CPS::GeneratorType::PVNode);
	systemPF.component<CPS::SP::Ph1::SynchronGenerator>("GEN1")->modifyPowerFlowBusType(CPS::PowerflowBusType::VD);

	// define logging
    auto loggerPF = DPsim::DataLogger::make(simNamePF);
    for (auto node : systemPF.mNodes)
    {
        loggerPF->logAttribute(node->name() + ".V", node->attribute("v"));
    }

	// run powerflow
    Simulation simPF(simNamePF, Logger::Level::debug);
	simPF.setSystem(systemPF);
	simPF.setTimeStep(finalTime);
	simPF.setFinalTime(2*finalTime);
	simPF.setDomain(Domain::SP);
	simPF.setSolverType(Solver::Type::NRP);
	simPF.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
	simPF.doInitFromNodesAndTerminals(true);
    simPF.addLogger(loggerPF);
    simPF.run();

	// ----- DYNAMIC SIMULATION -----
	Logger::setLogDir("logs/"+simName);

	CPS::CIM::Reader reader2(simName, Logger::Level::debug, Logger::Level::debug);
	SystemTopology sys = reader2.loadCIM(60, filenames, Domain::EMT, PhaseType::ABC, CPS::GeneratorType::FullOrder);

	sys.initWithPowerflow(systemPF);
	for (auto comp : sys.mComponents) {
		if (auto genEMT = std::dynamic_pointer_cast<CPS::EMT::Ph3::SynchronGeneratorDQTrapez>(comp)) {
			auto genPF = systemPF.component<CPS::SP::Ph1::SynchronGenerator>(comp->name());
			genEMT->terminal(0)->setPower(-genPF->getApparentPower());
		}
	}

	// Logging
	auto logger = DataLogger::make(simName);
	logger->logAttribute("v1", sys.node<SimNode>("BUS1")->attribute("v"));
	logger->logAttribute("v2", sys.node<SimNode>("BUS2")->attribute("v"));
	logger->logAttribute("v3", sys.node<SimNode>("BUS3")->attribute("v"));
	logger->logAttribute("v4", sys.node<SimNode>("BUS4")->attribute("v"));
	logger->logAttribute("v5", sys.node<SimNode>("BUS5")->attribute("v"));
	logger->logAttribute("v6", sys.node<SimNode>("BUS6")->attribute("v"));
	logger->logAttribute("v7", sys.node<SimNode>("BUS7")->attribute("v"));
	logger->logAttribute("v8", sys.node<SimNode>("BUS8")->attribute("v"));
	logger->logAttribute("v9", sys.node<SimNode>("BUS9")->attribute("v"));

	// log generator's current
	for (auto comp : sys.mComponents) {
		if (std::dynamic_pointer_cast<CPS::EMT::Ph3::SynchronGeneratorDQTrapez>(comp))
			logger->logAttribute(comp->name() + ".I", comp->attribute("i_intf"));
	}

	Simulation sim(simName, Logger::Level::info);
	sim.setSystem(sys);
	sim.setDomain(Domain::EMT);
	sim.setSolverType(Solver::Type::MNA);
	sim.setDirectLinearSolverImplementation(DPsim::DirectLinearSolverImpl::SparseLU);
	sim.setTimeStep(timeStep);
	sim.setFinalTime(finalTime);
	sim.addLogger(logger);
	sim.run();

	return 0;
}
