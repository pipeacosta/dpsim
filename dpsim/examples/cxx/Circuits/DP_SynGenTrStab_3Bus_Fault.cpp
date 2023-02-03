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
using namespace CIM::Examples::Grids::ThreeBus;

ScenarioConfig ThreeBus;

//Switch to trigger fault at generator terminal
Real SwitchOpen = 1e12;
Real SwitchClosed = 0.1;

void DP_SynGenTrStab_3Bus_Fault(String simName, Real timeStep, Real finalTime, bool startFaultEvent, bool endFaultEvent, Real startTimeFault, Real endTimeFault, Real cmdInertia_G1, Real cmdInertia_G2, Real cmdDamping_G1, Real cmdDamping_G2) {
	// ----- POWERFLOW FOR INITIALIZATION -----
	Real timeStepPF = finalTime;
	Real finalTimePF = finalTime+timeStepPF;
	String simNamePF = simName + "_PF";
	Logger::setLogDir("logs/" + simNamePF);

	// Components
	auto n1PF = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2PF = SimNode<Complex>::make("n2", PhaseType::Single);
	auto n3PF = SimNode<Complex>::make("n3", PhaseType::Single);

	//Synchronous generator 1
	auto gen1PF = SP::Ph1::SynchronGenerator::make("Generator", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	gen1PF->setParameters(ThreeBus.nomPower_G1, ThreeBus.nomPhPhVoltRMS_G1, ThreeBus.initActivePower_G1, ThreeBus.setPointVoltage_G1*ThreeBus.t1_ratio, PowerflowBusType::VD);
	gen1PF->setBaseVoltage(ThreeBus.Vnom);

	//Synchronous generator 2
	auto gen2PF = SP::Ph1::SynchronGenerator::make("Generator", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	gen2PF->setParameters(ThreeBus.nomPower_G2, ThreeBus.nomPhPhVoltRMS_G2, ThreeBus.initActivePower_G2, ThreeBus.setPointVoltage_G2*ThreeBus.t2_ratio, PowerflowBusType::PV);
	gen2PF->setBaseVoltage(ThreeBus.Vnom);

	//use Shunt as Load for powerflow
	auto loadPF = SP::Ph1::Shunt::make("Load", Logger::Level::debug);
	loadPF->setParameters(ThreeBus.activePower_L / std::pow(ThreeBus.Vnom, 2), - ThreeBus.reactivePower_L / std::pow(ThreeBus.Vnom, 2));
	loadPF->setBaseVoltage(ThreeBus.Vnom);

	//Line12
	auto line12PF = SP::Ph1::PiLine::make("PiLine12", Logger::Level::debug);
	line12PF->setParameters(ThreeBus.lineResistance12, ThreeBus.lineInductance12, ThreeBus.lineCapacitance12, ThreeBus.lineConductance12);
	line12PF->setBaseVoltage(ThreeBus.Vnom);
	//Line13
	auto line13PF = SP::Ph1::PiLine::make("PiLine13", Logger::Level::debug);
	line13PF->setParameters(ThreeBus.lineResistance13, ThreeBus.lineInductance13, ThreeBus.lineCapacitance13, ThreeBus.lineConductance13);
	line13PF->setBaseVoltage(ThreeBus.Vnom);
	//Line23
	auto line23PF = SP::Ph1::PiLine::make("PiLine23", Logger::Level::debug);
	line23PF->setParameters(ThreeBus.lineResistance23, ThreeBus.lineInductance23, ThreeBus.lineCapacitance23, ThreeBus.lineConductance23);
	line23PF->setBaseVoltage(ThreeBus.Vnom);
	//Switch
	auto faultPF = CPS::SP::Ph1::Switch::make("Br_fault", Logger::Level::debug);
	faultPF->setParameters(SwitchOpen, SwitchClosed);
	faultPF->open();

	// Topology
	gen1PF->connect({ n1PF });
	gen2PF->connect({ n2PF });
	loadPF->connect({ n3PF });
	line12PF->connect({ n1PF, n2PF });
	line13PF->connect({ n1PF, n3PF });
	line23PF->connect({ n2PF, n3PF });
	// faultPF->connect({SP::SimNode::GND , n1PF }); //terminal of generator 1
	faultPF->connect({SP::SimNode::GND , n2PF }); //terminal of generator 2
	// faultPF->connect({SP::SimNode::GND , n3PF }); //Load bus
	auto systemPF = SystemTopology(60,
			SystemNodeList{n1PF, n2PF, n3PF},
			SystemComponentList{gen1PF, gen2PF, loadPF, line12PF, line13PF, line23PF, faultPF});

	// Logging
	auto loggerPF = DataLogger::make(simNamePF);
	loggerPF->logAttribute("v_bus1", n1PF->attribute("v"));
	loggerPF->logAttribute("v_bus2", n2PF->attribute("v"));
	loggerPF->logAttribute("v_bus3", n3PF->attribute("v"));

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
	String simNameDP = simName + "_DP";
	Logger::setLogDir("logs/"+simNameDP);

	// Nodes
	auto n1DP = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2DP = SimNode<Complex>::make("n2", PhaseType::Single);
	auto n3DP = SimNode<Complex>::make("n3", PhaseType::Single);

	// Components
	//Synchronous generator 1
	auto gen1DP = DP::Ph1::SynchronGeneratorTrStab::make("SynGen1", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	gen1DP->setStandardParametersPU(ThreeBus.nomPower_G1, ThreeBus.nomPhPhVoltRMS_G1, ThreeBus.nomFreq_G1, ThreeBus.Xpd_G1*std::pow(ThreeBus.t1_ratio,2), cmdInertia_G1*ThreeBus.H_G1, ThreeBus.Rs_G1, cmdDamping_G1*ThreeBus.D_G1);
	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower_G1= gen1PF->getApparentPower();
	gen1DP->setInitialValues(initApparentPower_G1, ThreeBus.initMechPower_G1);

	//Synchronous generator 2
	auto gen2DP = DP::Ph1::SynchronGeneratorTrStab::make("SynGen2", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	gen2DP->setStandardParametersPU(ThreeBus.nomPower_G2, ThreeBus.nomPhPhVoltRMS_G2, ThreeBus.nomFreq_G2, ThreeBus.Xpd_G2*std::pow(ThreeBus.t2_ratio,2), cmdInertia_G2*ThreeBus.H_G2, ThreeBus.Rs_G2, cmdDamping_G2*ThreeBus.D_G2);
	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower_G2= gen2PF->getApparentPower();
	gen2DP->setInitialValues(initApparentPower_G2, ThreeBus.initMechPower_G2);

	gen2DP->setModelFlags(true);
	gen2DP->setReferenceOmega(gen1DP->attributeTyped<Real>("w_r"), gen1DP->attributeTyped<Real>("delta_r"));

	///Load
	auto loadDP=DP::Ph1::RXLoad::make("Load", Logger::Level::debug);
	loadDP->setParameters(ThreeBus.activePower_L, ThreeBus.reactivePower_L, ThreeBus.Vnom);

	//Line12
	auto line12DP = DP::Ph1::PiLine::make("PiLine12", Logger::Level::debug);
	line12DP->setParameters(ThreeBus.lineResistance12, ThreeBus.lineInductance12, ThreeBus.lineCapacitance12, ThreeBus.lineConductance12);
	//Line13
	auto line13DP = DP::Ph1::PiLine::make("PiLine13", Logger::Level::debug);
	line13DP->setParameters(ThreeBus.lineResistance13, ThreeBus.lineInductance13, ThreeBus.lineCapacitance13, ThreeBus.lineConductance13);
	//Line23
	auto line23DP = DP::Ph1::PiLine::make("PiLine23", Logger::Level::debug);
	line23DP->setParameters(ThreeBus.lineResistance23, ThreeBus.lineInductance23, ThreeBus.lineCapacitance23, ThreeBus.lineConductance23);

	// //Switch
	// auto faultDP = DP::Ph1::Switch::make("Br_fault", Logger::Level::debug);
	// faultDP->setParameters(SwitchOpen, SwitchClosed);
	// faultDP->open();

	//Variable resistance Switch
	auto faultDP = DP::Ph1::varResSwitch::make("Br_fault", Logger::Level::debug);
	faultDP->setParameters(SwitchOpen, SwitchClosed);
	faultDP->setInitParameters(timeStep);
	faultDP->open();

	// Topology
	gen1DP->connect({ n1DP });
	gen2DP->connect({ n2DP });
	loadDP->connect({ n3DP });
	line12DP->connect({ n1DP, n2DP });
	line13DP->connect({ n1DP, n3DP });
	line23DP->connect({ n2DP, n3DP });
	// faultDP->connect({DP::SimNode::GND , n1DP }); //terminal of generator 1
	faultDP->connect({DP::SimNode::GND , n2DP }); //terminal of generator 2
	// faultDP->connect({DP::SimNode::GND , n3DP }); //Load bus
	auto systemDP = SystemTopology(60,
			SystemNodeList{n1DP, n2DP, n3DP},
			SystemComponentList{gen1DP, gen2DP, loadDP, line12DP, line13DP, line23DP, faultDP});

	// Initialization of dynamic topology
	systemDP.initWithPowerflow(systemPF);

	// Logging
	auto loggerDP = DataLogger::make(simNameDP);
	loggerDP->logAttribute("v1", n1DP->attribute("v"));
	loggerDP->logAttribute("v2", n2DP->attribute("v"));
	loggerDP->logAttribute("v3", n3DP->attribute("v"));
	loggerDP->logAttribute("v_line12", line12DP->attribute("v_intf"));
	loggerDP->logAttribute("i_line12", line12DP->attribute("i_intf"));
	loggerDP->logAttribute("v_line13", line13DP->attribute("v_intf"));
	loggerDP->logAttribute("i_line13", line13DP->attribute("i_intf"));
	loggerDP->logAttribute("v_line23", line23DP->attribute("v_intf"));
	loggerDP->logAttribute("i_line23", line23DP->attribute("i_intf"));
	loggerDP->logAttribute("Ep_gen1", gen1DP->attribute("Ep_mag"));
	loggerDP->logAttribute("v_gen1", gen1DP->attribute("v_intf"));
	loggerDP->logAttribute("i_gen1", gen1DP->attribute("i_intf"));
	loggerDP->logAttribute("wr_gen1", gen1DP->attribute("w_r"));
	loggerDP->logAttribute("delta_gen1", gen1DP->attribute("delta_r"));
	loggerDP->logAttribute("Ep_gen2", gen2DP->attribute("Ep_mag"));
	loggerDP->logAttribute("v_gen2", gen2DP->attribute("v_intf"));
	loggerDP->logAttribute("i_gen2", gen2DP->attribute("i_intf"));
	loggerDP->logAttribute("wr_gen2", gen2DP->attribute("w_r"));
	loggerDP->logAttribute("wref_gen2", gen2DP->attribute("w_ref"));
	loggerDP->logAttribute("delta_gen2", gen2DP->attribute("delta_r"));
	loggerDP->logAttribute("i_fault", faultDP->attribute("i_intf"));
	loggerDP->logAttribute("v_load", loadDP->attribute("v_intf"));
	loggerDP->logAttribute("i_load", loadDP->attribute("i_intf"));
	loggerDP->logAttribute("P_mech1", gen1DP->attribute("P_mech"));
	loggerDP->logAttribute("P_mech2", gen2DP->attribute("P_mech"));
	loggerDP->logAttribute("P_elec1", gen1DP->attribute("P_elec"));
	loggerDP->logAttribute("P_elec2", gen2DP->attribute("P_elec"));

	Simulation simDP(simNameDP, Logger::Level::debug);
	simDP.setSystem(systemDP);
	simDP.setTimeStep(timeStep);
	simDP.setFinalTime(finalTime);
	simDP.setDomain(Domain::DP);
	simDP.addLogger(loggerDP);
	simDP.doSystemMatrixRecomputation(true);
	simDP.setMnaSolverImplementation(MnaSolverFactory::MnaSolverImpl::EigenSparse);

	// Events
	if (startFaultEvent){
		auto sw1 = SwitchEvent::make(startTimeFault, faultDP, true);

		simDP.addEvent(sw1);
	}

	if(endFaultEvent){

		auto sw2 = SwitchEvent::make(endTimeFault, faultDP, false);
		simDP.addEvent(sw2);

	}

	simDP.run();
}

int main(int argc, char* argv[]) {


	//Simultion parameters
	String simName="DP_SynGenTrStab_3Bus_Fault";
	Real finalTime = 30;
	Real timeStep = 0.001;
	Bool startFaultEvent=true;
	Bool endFaultEvent=true;
	Real startTimeFault=10;
	Real endTimeFault=10.2;
	Real cmdInertia_G1= 1.0;
	Real cmdInertia_G2= 1.0;
	Real cmdDamping_G1=1.0;
	Real cmdDamping_G2=1.0;

	CommandLineArgs args(argc, argv);
	if (argc > 1) {
		timeStep = args.timeStep;
		finalTime = args.duration;
		if (args.name != "dpsim")
			simName = args.name;
		if (args.options.find("SCALEINERTIA_G1") != args.options.end())
			cmdInertia_G1 = args.getOptionReal("SCALEINERTIA_G1");
		if (args.options.find("SCALEINERTIA_G2") != args.options.end())
			cmdInertia_G2 = args.getOptionReal("SCALEINERTIA_G2");
		if (args.options.find("SCALEDAMPING_G1") != args.options.end())
			cmdDamping_G1 = args.getOptionReal("SCALEDAMPING_G1");
		if (args.options.find("SCALEDAMPING_G2") != args.options.end())
			cmdDamping_G2 = args.getOptionReal("SCALEDAMPING_G2");
		if (args.options.find("STARTTIMEFAULT") != args.options.end())
			startTimeFault = args.getOptionReal("STARTTIMEFAULT");
		if (args.options.find("ENDTIMEFAULT") != args.options.end())
			endTimeFault = args.getOptionReal("ENDTIMEFAULT");
	}

	DP_SynGenTrStab_3Bus_Fault(simName, timeStep, finalTime, startFaultEvent, endFaultEvent, startTimeFault, endTimeFault, cmdInertia_G1, cmdInertia_G2, cmdDamping_G1, cmdDamping_G2);
}
