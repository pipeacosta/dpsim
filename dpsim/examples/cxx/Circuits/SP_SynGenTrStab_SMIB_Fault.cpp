/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <DPsim.h>
#include "../Examples.h"

using namespace DPsim;
using namespace CPS;
using namespace CIM::Examples::Grids::SMIB;

ScenarioConfig smib;

//Switch to trigger fault at generator terminal
Real SwitchOpen = 1e6;
Real SwitchClosed = 0.1;

void SP_1ph_SynGenTrStab_Fault(String simName, Real timeStep, Real finalTime, bool startFaultEvent, bool endFaultEvent, Real startTimeFault, Real endTimeFault, Real cmdInertia, Real cmdDamping) {
	// ----- POWERFLOW FOR INITIALIZATION -----
	Real timeStepPF = finalTime;
	Real finalTimePF = finalTime+timeStepPF;
	String simNamePF = simName + "_PF";
	Logger::setLogDir("logs/" + simNamePF);

	// Components
	auto n1PF = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2PF = SimNode<Complex>::make("n2", PhaseType::Single);

	//Synchronous generator ideal model
	auto genPF = SP::Ph1::SynchronGenerator::make("Generator", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	genPF->setParameters(smib.nomPower, smib.nomPhPhVoltRMS, smib.initActivePower, smib.setPointVoltage*smib.t_ratio, PowerflowBusType::PV);
	genPF->setBaseVoltage(smib.Vnom);
	genPF->modifyPowerFlowBusType(PowerflowBusType::PV);

	//Grid bus as Slack
	auto extnetPF = SP::Ph1::NetworkInjection::make("Slack", Logger::Level::debug);
	extnetPF->setParameters(smib.Vnom);
	extnetPF->setBaseVoltage(smib.Vnom);
	extnetPF->modifyPowerFlowBusType(PowerflowBusType::VD);

	//Line
	auto linePF = SP::Ph1::PiLine::make("PiLine", Logger::Level::debug);
	linePF->setParameters(smib.lineResistance, smib.lineInductance, smib.lineCapacitance, smib.lineConductance);
	linePF->setBaseVoltage(smib.Vnom);

	//Switch
	auto faultPF = CPS::SP::Ph1::Switch::make("Br_fault", Logger::Level::debug);
	faultPF->setParameters(SwitchOpen, SwitchClosed);
	faultPF->open();

	// Topology
	genPF->connect({ n1PF });
	faultPF->connect({SP::SimNode::GND , n1PF });
	linePF->connect({ n1PF, n2PF });
	extnetPF->connect({ n2PF });
	auto systemPF = SystemTopology(60,
			SystemNodeList{n1PF, n2PF},
			SystemComponentList{genPF, linePF, extnetPF, faultPF});

	// Logging
	auto loggerPF = DataLogger::make(simNamePF);
	loggerPF->logAttribute("v1", n1PF->attribute("v"));
	loggerPF->logAttribute("v2", n2PF->attribute("v"));

	// Simulation
	Simulation simPF(simNamePF, Logger::Level::debug);
	simPF.setSystem(systemPF);
	simPF.setTimeStep(timeStepPF);
	simPF.setFinalTime(finalTimePF);
	simPF.setDomain(Domain::SP);
	simPF.setSolverType(Solver::Type::NRP);
	simPF.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
	simPF.doInitFromNodesAndTerminals(false);
	simPF.addLogger(loggerPF);
	simPF.run();

	// ----- Dynamic simulation ------
	String simNameSP = simName + "_SP";
	Logger::setLogDir("logs/"+simNameSP);

	// Nodes
	auto n1SP = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2SP = SimNode<Complex>::make("n2", PhaseType::Single);

	// Components
	auto genSP = CPS::SP::Ph1::SynchronGeneratorTrStab::make("SynGen", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	genSP->setStandardParametersPU(smib.nomPower, smib.nomPhPhVoltRMS, smib.nomFreq, smib.Xpd*std::pow(smib.t_ratio,2), cmdInertia*smib.H, smib.Rs, cmdDamping*smib.D );
	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower= genPF->getApparentPower();
	genSP->setInitialValues(initApparentPower, smib.initMechPower);

	//Grid bus as Slack
	auto extnetSP = SP::Ph1::NetworkInjection::make("Slack", Logger::Level::debug);
	extnetSP->setParameters(smib.Vnom);
	// Line
	auto lineSP = SP::Ph1::PiLine::make("PiLine", Logger::Level::debug);
	lineSP->setParameters(smib.lineResistance, smib.lineInductance, smib.lineCapacitance, smib.lineConductance);

	// //Switch
	// auto faultSP = SP::Ph1::Switch::make("Br_fault", Logger::Level::debug);
	// faultSP->setParameters(SwitchOpen, SwitchClosed);
	// faultSP->open();

	// Variable resistance Switch
	auto faultSP = SP::Ph1::varResSwitch::make("Br_fault", Logger::Level::debug);
	faultSP->setParameters(SwitchOpen, SwitchClosed);
	faultSP->setInitParameters(timeStep);
	faultSP->open();

	// Topology
	genSP->connect({ n1SP });
	faultSP->connect({SP::SimNode::GND , n1SP });
	lineSP->connect({ n1SP, n2SP });
	extnetSP->connect({ n2SP });
	auto systemSP = SystemTopology(60,
			SystemNodeList{n1SP, n2SP},
			SystemComponentList{genSP, lineSP, extnetSP, faultSP});

	// Initialization of dynamic topology
	systemSP.initWithPowerflow(systemPF);


	// Logging
	auto loggerSP = DataLogger::make(simNameSP);
	loggerSP->logAttribute("v1", n1SP->attribute("v"));
	loggerSP->logAttribute("v2", n2SP->attribute("v"));
	//gen
	loggerSP->logAttribute("Ep", genSP->attribute("Ep"));
	loggerSP->logAttribute("v_gen", genSP->attribute("v_intf"));
	loggerSP->logAttribute("i_gen", genSP->attribute("i_intf"));
	loggerSP->logAttribute("wr_gen", genSP->attribute("w_r"));
	loggerSP->logAttribute("delta_r_gen", genSP->attribute("delta_r"));
	loggerSP->logAttribute("P_elec", genSP->attribute("P_elec"));
	loggerSP->logAttribute("P_mech", genSP->attribute("P_mech"));
	//Switch
	loggerSP->logAttribute("i_fault", faultSP->attribute("i_intf"));
	//line
	loggerSP->logAttribute("v_line", lineSP->attribute("v_intf"));
	loggerSP->logAttribute("i_line", lineSP->attribute("i_intf"));
	//slack
	loggerSP->logAttribute("v_slack", extnetSP->attribute("v_intf"));
	loggerSP->logAttribute("i_slack", extnetSP->attribute("i_intf"));


	Simulation simSP(simNameSP, Logger::Level::debug);
	simSP.setSystem(systemSP);
	simSP.setTimeStep(timeStep);
	simSP.setFinalTime(finalTime);
	simSP.setDomain(Domain::SP);
	simSP.addLogger(loggerSP);
	simSP.doSystemMatrixRecomputation(true);
	simSP.setDirectLinearSolverImplementation(DPsim::DirectLinearSolverImpl::SparseLU);

	// Events
	if (startFaultEvent){
		auto sw1 = SwitchEvent::make(startTimeFault, faultSP, true);

		simSP.addEvent(sw1);
	}

	if(endFaultEvent){

		auto sw2 = SwitchEvent::make(endTimeFault, faultSP, false);
		simSP.addEvent(sw2);

	}

	simSP.run();
}

int main(int argc, char* argv[]) {


	//Simultion parameters
	String simName="SP_SynGenTrStab_SMIB_Fault";
	Real finalTime = 30;
	Real timeStep = 0.001;
	Bool startFaultEvent=true;
	Bool endFaultEvent=true;
	Real startTimeFault=10;
	Real endTimeFault=10.2;
	Real cmdInertia= 1.0;
	Real cmdDamping=1.0;

	CommandLineArgs args(argc, argv);
	if (argc > 1) {
		timeStep = args.timeStep;
		finalTime = args.duration;
		if (args.name != "dpsim")
			simName = args.name;
		if (args.options.find("SCALEINERTIA") != args.options.end())
			cmdInertia = args.getOptionReal("SCALEINERTIA");
		if (args.options.find("SCALEDAMPING") != args.options.end())
			cmdDamping = args.getOptionReal("SCALEDAMPING");
		if (args.options.find("STARTTIMEFAULT") != args.options.end())
			startTimeFault = args.getOptionReal("STARTTIMEFAULT");
		if (args.options.find("ENDTIMEFAULT") != args.options.end())
			endTimeFault = args.getOptionReal("ENDTIMEFAULT");
	}

	SP_1ph_SynGenTrStab_Fault(simName, timeStep, finalTime, startFaultEvent, endFaultEvent, startTimeFault, endTimeFault, cmdInertia, cmdDamping);
}
