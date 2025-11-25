/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University; Universidad
 *                     Nacional de Colombia
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <DPsim.h>
#include "../Examples.h"
#include "dpsim-models/Definitions.h"
#include "dpsim-models/IdentifiedObject.h"
#include "dpsim-models/SP/SP_Ph1_PiLine.h"
#include "dpsim-models/Signal/DecouplingIdealTransformer_SP_Ph1.h"
#include "dpsim-models/SimNode.h"
#include "dpsim/Definitions.h"

using namespace DPsim;
using namespace CPS;
using namespace CIM::Examples::Grids::KRK_TwoArea;
using namespace CIM::Examples;

ScenarioConfig KRK_TwoArea;

void decoupleNode(SystemTopology &sys, const String &nodeName, const IdentifiedObject::List &componentsAt1,
                    const IdentifiedObject::List &componentsAt2, Real delay, CouplingMethod method) {
  SimNode<Complex>::List newNodes;
  SimPowerComp<Complex>::List newComponents;

  auto intfNode = sys.node<SP::SimNode>(nodeName);
  std::shared_ptr<TopologicalNode> nodeCopy1Topo = intfNode->clone(nodeName + "_1");
  std::shared_ptr<TopologicalNode> nodeCopy2Topo = intfNode->clone(nodeName + "_2");

  auto nodeCopy1 = std::dynamic_pointer_cast<SimNode<Complex>>(nodeCopy1Topo);
  auto nodeCopy2 = std::dynamic_pointer_cast<SimNode<Complex>>(nodeCopy2Topo);

  newNodes.push_back(nodeCopy1);
  newNodes.push_back(nodeCopy2);

  for (auto genComp : componentsAt1) {
    auto comp = std::dynamic_pointer_cast<SimPowerComp<Complex>>(genComp);
    auto compCopy = comp->clone(comp->name());

    SimNode<Complex>::List nodeCopies;
    for (UInt nNode = 0; nNode < comp->terminalNumber(); nNode++) {
      String origNodeName = comp->node(nNode)->name();
      if (origNodeName == nodeName) {
        nodeCopies.push_back(nodeCopy1);
      } else {
        nodeCopies.push_back(comp->node(nNode));
      }
    }

    compCopy->connect(nodeCopies);

    // update the terminal powers for powerflow initialization
    for (UInt nTerminal = 0; nTerminal < comp->terminalNumber();
          nTerminal++) {
      compCopy->terminal(nTerminal)->setPower(comp->terminal(nTerminal)->power());
    }
    newComponents.push_back(compCopy);
    sys.removeComponent(comp->name());
  }

  for (auto genComp : componentsAt2) {
    auto comp = std::dynamic_pointer_cast<SimPowerComp<Complex>>(genComp);
    auto compCopy = comp->clone(comp->name());

    SimNode<Complex>::List nodeCopies;
    for (UInt nNode = 0; nNode < comp->terminalNumber(); nNode++) {
      String origNodeName = comp->node(nNode)->name();
      if (origNodeName == nodeName) {
        nodeCopies.push_back(nodeCopy2);
      } else {
        nodeCopies.push_back(comp->node(nNode));
      }
    }

    compCopy->connect(nodeCopies);

    // update the terminal powers for powerflow initialization
    for (UInt nTerminal = 0; nTerminal < comp->terminalNumber();
          nTerminal++) {
      compCopy->terminal(nTerminal)->setPower(comp->terminal(nTerminal)->power());
    }
    newComponents.push_back(compCopy);
    sys.removeComponent(comp->name());
  }

  sys.removeNode(nodeName);

  for (auto node : newNodes)
    sys.addNode(node);
  for (auto comp : newComponents) {
    sys.addComponent(comp);
  }

  Matrix i_0 = Matrix::Zero(3,1);

  auto idealTrafo = Signal::DecouplingIdealTransformer_SP_Ph1::make("itm_" + nodeName,
                                                                Logger::Level::debug);
  idealTrafo->setParameters(nodeCopy1, nodeCopy2, delay, i_0, Complex(0, 0), method);
  sys.addComponent(idealTrafo);
  sys.addComponents(idealTrafo->getComponents());
  sys.addNode(idealTrafo->getVirtualNode());
}

void SP_SynGenTrStab_KRK_TwoAreaTrafo_SteadyState(String simName, Real timeStep, Real finalTime, bool startFaultEvent, bool endFaultEvent, Real startTimeFault, Real endTimeFault, Real cmdInertia_G1, Real cmdInertia_G2, Real cmdInertia_G3, Real cmdInertia_G4, Real cmdDamping_G1, Real cmdDamping_G2, Real cmdDamping_G3, Real cmdDamping_G4, Real delay, CouplingMethod method) {
	// ----- POWERFLOW FOR INITIALIZATION -----
	Real timeStepPF = finalTime;
	Real finalTimePF = finalTime+ timeStepPF;
	String simNamePF = simName + "_PF";
	Logger::setLogDir("logs/" + simNamePF);

	// Components
	auto n1PF = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2PF = SimNode<Complex>::make("n2", PhaseType::Single);
	auto n3PF = SimNode<Complex>::make("n3", PhaseType::Single);
  auto n4PF = SimNode<Complex>::make("n4", PhaseType::Single);
	auto n5PF = SimNode<Complex>::make("n5", PhaseType::Single);
	auto n6PF = SimNode<Complex>::make("n6", PhaseType::Single);
  auto n7PF = SimNode<Complex>::make("n7", PhaseType::Single);
	auto n8PF = SimNode<Complex>::make("n8", PhaseType::Single);
	auto n9PF = SimNode<Complex>::make("n9", PhaseType::Single);
  auto n10PF = SimNode<Complex>::make("n10", PhaseType::Single);
	auto n11PF = SimNode<Complex>::make("n11", PhaseType::Single);

	//Synchronous generator 1
	auto gen1PF = SP::Ph1::SynchronGenerator::make("SynGen1", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	gen1PF->setParameters(KRK_TwoArea.nomPower_G1, KRK_TwoArea.nomPhPhVoltRMS_G1, KRK_TwoArea.initActivePower_G1, KRK_TwoArea.setPointVoltage_G1, PowerflowBusType::PV);
	gen1PF->setBaseVoltage(KRK_TwoArea.nomPhPhVoltRMS_G1);

	//Synchronous generator 2
	auto gen2PF = SP::Ph1::SynchronGenerator::make("SynGen2", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	gen2PF->setParameters(KRK_TwoArea.nomPower_G2, KRK_TwoArea.nomPhPhVoltRMS_G2, KRK_TwoArea.initActivePower_G2, KRK_TwoArea.setPointVoltage_G2, PowerflowBusType::PV);
	gen2PF->setBaseVoltage(KRK_TwoArea.nomPhPhVoltRMS_G2);

    //Synchronous generator 3
	auto gen3PF = SP::Ph1::SynchronGenerator::make("SynGen3", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	gen3PF->setParameters(KRK_TwoArea.nomPower_G3, KRK_TwoArea.nomPhPhVoltRMS_G3, KRK_TwoArea.initActivePower_G3, KRK_TwoArea.setPointVoltage_G3, PowerflowBusType::VD);
	gen3PF->setBaseVoltage(KRK_TwoArea.nomPhPhVoltRMS_G3);

    //Synchronous generator 4
	auto gen4PF = SP::Ph1::SynchronGenerator::make("SynGen4", Logger::Level::debug);
	// setPointVoltage is defined as the voltage at the transfomer primary side and should be transformed to network side
	gen4PF->setParameters(KRK_TwoArea.nomPower_G4, KRK_TwoArea.nomPhPhVoltRMS_G4, KRK_TwoArea.initActivePower_G4, KRK_TwoArea.setPointVoltage_G4, PowerflowBusType::PV);
	gen4PF->setBaseVoltage(KRK_TwoArea.nomPhPhVoltRMS_G4);

	auto trafo15PF= SP::Ph1::Transformer::make("Transformer15", Logger::Level::debug);
	auto trafo26PF= SP::Ph1::Transformer::make("Transformer26", Logger::Level::debug);
	auto trafo311PF= SP::Ph1::Transformer::make("Transformer311", Logger::Level::debug);
	auto trafo410PF= SP::Ph1::Transformer::make("Transformer410", Logger::Level::debug);

	Real voltageMVSide = KRK_TwoArea.nomPhPhVoltRMS_G1;
	Real voltageHVSide = KRK_TwoArea.Vnom;
	Real ratio = voltageMVSide / voltageHVSide;
	Real trafoResistance = 0;
	Real trafoInductance = 23.4e-3; //set base to V_end2
	Real trafoPower = 900e6;
	// Note: to be consistent impedance values must be referred to high voltage side (and base voltage set to higher voltage)
	trafo15PF->setParameters(voltageMVSide, voltageHVSide, trafoPower, std::abs(ratio), std::arg(ratio), trafoResistance, trafoInductance);
	trafo15PF->setBaseVoltage(KRK_TwoArea.Vnom);

	trafo26PF->setParameters(voltageMVSide, voltageHVSide, trafoPower, std::abs(ratio), std::arg(ratio), trafoResistance, trafoInductance);
	trafo26PF->setBaseVoltage(KRK_TwoArea.Vnom);

	trafo311PF->setParameters(voltageMVSide, voltageHVSide, trafoPower, std::abs(ratio), std::arg(ratio), trafoResistance, trafoInductance);
	trafo311PF->setBaseVoltage(KRK_TwoArea.Vnom);

	trafo410PF->setParameters(voltageMVSide, voltageHVSide, trafoPower, std::abs(ratio), std::arg(ratio), trafoResistance, trafoInductance);
	trafo410PF->setBaseVoltage(KRK_TwoArea.Vnom);

    //use Shunt as Load for powerflow
	auto load7PF = SP::Ph1::Load::make("Load7", Logger::Level::debug);
	load7PF->setParameters(KRK_TwoArea.activePower_L7, KRK_TwoArea.reactivePower_L7_inductive - KRK_TwoArea.reactivePower_L7_capacitive, KRK_TwoArea.Vnom);

	auto load9PF = SP::Ph1::Load::make("Load9", Logger::Level::debug);
	load9PF->setParameters(KRK_TwoArea.activePower_L9, KRK_TwoArea.reactivePower_L9_inductive - KRK_TwoArea.reactivePower_L9_capacitive, KRK_TwoArea.Vnom);

	//Line56
	auto line56PF = SP::Ph1::PiLine::make("PiLine56", Logger::Level::debug);
	line56PF->setParameters(KRK_TwoArea.lineResistance56, KRK_TwoArea.lineInductance56, KRK_TwoArea.lineCapacitance56, KRK_TwoArea.lineConductance56);
	line56PF->setBaseVoltage(KRK_TwoArea.Vnom);

  //Line67
	auto line67PF = SP::Ph1::PiLine::make("PiLine67", Logger::Level::debug);
	line67PF->setParameters(KRK_TwoArea.lineResistance67, KRK_TwoArea.lineInductance67, KRK_TwoArea.lineCapacitance67, KRK_TwoArea.lineConductance67);
	line67PF->setBaseVoltage(KRK_TwoArea.Vnom);

  //Line78_1
	auto line78_1PF = SP::Ph1::PiLine::make("Piline78_1", Logger::Level::debug);
	line78_1PF->setParameters(KRK_TwoArea.lineResistance78, KRK_TwoArea.lineInductance78, KRK_TwoArea.lineCapacitance78, KRK_TwoArea.lineConductance78);
	line78_1PF->setBaseVoltage(KRK_TwoArea.Vnom);

  //Line78_2
	auto line78_2PF = SP::Ph1::PiLine::make("Piline78_2", Logger::Level::debug);
	line78_2PF->setParameters(KRK_TwoArea.lineResistance78, KRK_TwoArea.lineInductance78, KRK_TwoArea.lineCapacitance78, KRK_TwoArea.lineConductance78);
	line78_2PF->setBaseVoltage(KRK_TwoArea.Vnom);

  //Line89_1
	auto line89_1PF = SP::Ph1::PiLine::make("Piline89_1", Logger::Level::debug);
	line89_1PF->setParameters(KRK_TwoArea.lineResistance89, KRK_TwoArea.lineInductance89, KRK_TwoArea.lineCapacitance89, KRK_TwoArea.lineConductance89);
	line89_1PF->setBaseVoltage(KRK_TwoArea.Vnom);

  //Line89_2
	auto line89_2PF = SP::Ph1::PiLine::make("Piline89_2", Logger::Level::debug);
	line89_2PF->setParameters(KRK_TwoArea.lineResistance89, KRK_TwoArea.lineInductance89, KRK_TwoArea.lineCapacitance89, KRK_TwoArea.lineConductance89);
	line89_2PF->setBaseVoltage(KRK_TwoArea.Vnom);

  //Line910
	auto line910PF = SP::Ph1::PiLine::make("PiLine910", Logger::Level::debug);
	line910PF->setParameters(KRK_TwoArea.lineResistance910, KRK_TwoArea.lineInductance910, KRK_TwoArea.lineCapacitance910, KRK_TwoArea.lineConductance910);
	line910PF->setBaseVoltage(KRK_TwoArea.Vnom);

  //Line1011
	auto line1011PF = SP::Ph1::PiLine::make("PiLine1011", Logger::Level::debug);
	line1011PF->setParameters(KRK_TwoArea.lineResistance1011, KRK_TwoArea.lineInductance1011, KRK_TwoArea.lineCapacitance1011, KRK_TwoArea.lineConductance1011);
	line1011PF->setBaseVoltage(KRK_TwoArea.Vnom);

	// Topology
	gen1PF->connect({ n1PF });
	gen2PF->connect({ n2PF });
  gen3PF->connect({ n3PF });
  gen4PF->connect({ n4PF });

	trafo15PF->connect({ n1PF, n5PF});
	trafo26PF->connect({ n2PF, n6PF});
	trafo311PF->connect({ n3PF, n11PF});
	trafo410PF->connect({ n4PF, n10PF});

	load7PF->connect({ n7PF });
  load9PF->connect({ n9PF });

	line56PF->connect({ n5PF, n6PF });
	line67PF->connect({ n6PF, n7PF });
	line78_1PF->connect({ n7PF, n8PF });
	line78_2PF->connect({ n7PF, n8PF });
  line89_1PF->connect({ n8PF, n9PF });
	line89_2PF->connect({ n8PF, n9PF });
  line910PF->connect({ n9PF, n10PF });
  line1011PF->connect({ n10PF, n11PF });
	auto systemPF = SystemTopology(60,
			SystemNodeList{ n1PF, n2PF, n3PF, n4PF, n5PF, n6PF, n7PF, n8PF, n9PF, n10PF, n11PF},
			SystemComponentList{gen1PF, gen2PF, gen3PF, gen4PF, trafo15PF, trafo26PF, trafo311PF, trafo410PF, load7PF, load9PF, line56PF, line67PF, line78_1PF, line78_2PF, line89_1PF, line89_2PF, line910PF, line1011PF});

	// Logging
	auto loggerPF = DataLogger::make(simNamePF);
	loggerPF->logAttribute("v_bus1", n1PF->attribute("v"));
	loggerPF->logAttribute("s_bus1", n1PF->attribute("s"));
	loggerPF->logAttribute("v_bus2", n2PF->attribute("v"));
	loggerPF->logAttribute("s_bus2", n2PF->attribute("s"));
	loggerPF->logAttribute("v_bus3", n3PF->attribute("v"));
	loggerPF->logAttribute("s_bus3", n3PF->attribute("s"));
  loggerPF->logAttribute("v_bus4", n4PF->attribute("v"));
	loggerPF->logAttribute("s_bus4", n4PF->attribute("s"));
	loggerPF->logAttribute("v_bus5", n5PF->attribute("v"));
	loggerPF->logAttribute("s_bus5", n5PF->attribute("s"));
	loggerPF->logAttribute("v_bus6", n6PF->attribute("v"));
	loggerPF->logAttribute("s_bus6", n6PF->attribute("s"));
  loggerPF->logAttribute("v_bus7", n7PF->attribute("v"));
	loggerPF->logAttribute("s_bus7", n7PF->attribute("s"));
	loggerPF->logAttribute("v_bus8", n8PF->attribute("v"));
	loggerPF->logAttribute("s_bus8", n8PF->attribute("s"));
	loggerPF->logAttribute("v_bus9", n9PF->attribute("v"));
	loggerPF->logAttribute("s_bus9", n9PF->attribute("s"));
  loggerPF->logAttribute("v_bus10", n10PF->attribute("v"));
	loggerPF->logAttribute("s_bus10", n10PF->attribute("s"));
	loggerPF->logAttribute("v_bus11", n11PF->attribute("v"));
	loggerPF->logAttribute("s_bus11", n11PF->attribute("s"));

	// Simulation
	Simulation simPF(simNamePF, Logger::Level::debug);
	simPF.setSystem(systemPF);
	simPF.setTimeStep(timeStepPF);
	simPF.setFinalTime(finalTimePF);
	simPF.setDomain(Domain::SP);
	simPF.setSolverType(Solver::Type::NRP);
	simPF.doInitFromNodesAndTerminals(false);
	simPF.addLogger(loggerPF);
	simPF.run();

	// ----- Dynamic simulation ------
	String simNameSP = simName + "_SP";
	Logger::setLogDir("logs/"+simNameSP);

	// Nodes
	auto n1SP = SimNode<Complex>::make("n1", PhaseType::Single);
	auto n2SP = SimNode<Complex>::make("n2", PhaseType::Single);
	auto n3SP = SimNode<Complex>::make("n3", PhaseType::Single);
	auto n4SP = SimNode<Complex>::make("n4", PhaseType::Single);
	auto n5SP = SimNode<Complex>::make("n5", PhaseType::Single);
	auto n6SP = SimNode<Complex>::make("n6", PhaseType::Single);
	auto n7SP = SimNode<Complex>::make("n7", PhaseType::Single);
	auto n8SP = SimNode<Complex>::make("n8", PhaseType::Single);
	auto n9SP = SimNode<Complex>::make("n9", PhaseType::Single);
	auto n10SP = SimNode<Complex>::make("n10", PhaseType::Single);
	auto n11SP = SimNode<Complex>::make("n11", PhaseType::Single);

	// Components
	//Synchronous generator 1
	auto gen1SP = SP::Ph1::SynchronGeneratorTrStab::make("SynGen1", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	gen1SP->setStandardParametersPU(KRK_TwoArea.nomPower_G1, KRK_TwoArea.nomPhPhVoltRMS_G1, KRK_TwoArea.nomFreq_G1, KRK_TwoArea.Xpd, KRK_TwoArea.H_G1, KRK_TwoArea.Rs, 10.0);
	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower_G1 = gen1PF->getApparentPower();
	Real initMechPower_G1 = initApparentPower_G1.real();
	gen1SP->setInitialValues(initApparentPower_G1, initMechPower_G1);

	//Synchronous generator 2
	auto gen2SP = SP::Ph1::SynchronGeneratorTrStab::make("SynGen2", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	gen2SP->setStandardParametersPU(KRK_TwoArea.nomPower_G2, KRK_TwoArea.nomPhPhVoltRMS_G2, KRK_TwoArea.nomFreq_G2, KRK_TwoArea.Xpd, KRK_TwoArea.H_G2, KRK_TwoArea.Rs, 10.0);
	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower_G2 = gen2PF->getApparentPower();
	Real initMechPower_G2 = initApparentPower_G2.real();
	gen2SP->setInitialValues(initApparentPower_G2, initMechPower_G2);

	//Synchronous generator 3
	auto gen3SP = SP::Ph1::SynchronGeneratorTrStab::make("SynGen3", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	gen3SP->setStandardParametersPU(KRK_TwoArea.nomPower_G3, KRK_TwoArea.nomPhPhVoltRMS_G3, KRK_TwoArea.nomFreq_G3, KRK_TwoArea.Xpd, KRK_TwoArea.H_G3, KRK_TwoArea.Rs, 10.0);
	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower_G3 = gen3PF->getApparentPower();
	Real initMechPower_G3 = initApparentPower_G3.real();
	gen3SP->setInitialValues(initApparentPower_G3, initMechPower_G3);

	//Synchronous generator 4
	auto gen4SP = SP::Ph1::SynchronGeneratorTrStab::make("SynGen4", Logger::Level::debug);
	// Xpd is given in p.u of generator base at transfomer primary side and should be transformed to network side
	gen4SP->setStandardParametersPU(KRK_TwoArea.nomPower_G4, KRK_TwoArea.nomPhPhVoltRMS_G4, KRK_TwoArea.nomFreq_G4, KRK_TwoArea.Xpd, KRK_TwoArea.H_G4, KRK_TwoArea.Rs, 10.0);
	// Get actual active and reactive power of generator's Terminal from Powerflow solution
	Complex initApparentPower_G4 = gen4PF->getApparentPower();
	Real initMechPower_G4 = initApparentPower_G4.real();
	gen4SP->setInitialValues(initApparentPower_G4, initMechPower_G4);

	auto trafo15SP = SP::Ph1::Transformer::make("trafo15", Logger::Level::debug);
	trafo15SP->setParameters(voltageMVSide, voltageHVSide, trafoPower, std::abs(ratio), std::arg(ratio), trafoResistance, trafoInductance);

	auto trafo26SP = SP::Ph1::Transformer::make("trafo26", Logger::Level::debug);
	trafo26SP->setParameters(voltageMVSide, voltageHVSide, trafoPower, std::abs(ratio), std::arg(ratio), trafoResistance, trafoInductance);

	auto trafo311SP = SP::Ph1::Transformer::make("trafo311", Logger::Level::debug);
	trafo311SP->setParameters(voltageMVSide, voltageHVSide, trafoPower, std::abs(ratio), std::arg(ratio), trafoResistance, trafoInductance);

	auto trafo410SP = SP::Ph1::Transformer::make("trafo410", Logger::Level::debug);
	trafo410SP->setParameters(voltageMVSide, voltageHVSide, trafoPower, std::abs(ratio), std::arg(ratio), trafoResistance, trafoInductance);

	///Loads
	auto load7SP = SP::Ph1::Load::make("Load7", Logger::Level::debug);
	load7SP->setParameters(KRK_TwoArea.activePower_L7, KRK_TwoArea.reactivePower_L7_inductive - KRK_TwoArea.reactivePower_L7_capacitive, KRK_TwoArea.Vnom);

	auto load9SP = SP::Ph1::Load::make("Load9", Logger::Level::debug);
	load9SP->setParameters(KRK_TwoArea.activePower_L9, KRK_TwoArea.reactivePower_L9_inductive - KRK_TwoArea.reactivePower_L9_capacitive, KRK_TwoArea.Vnom);

	//Line45
	auto line45SP = SP::Ph1::PiLine::make("PiLine45", Logger::Level::debug);
	line45SP->setParameters(KRK_TwoArea.lineResistance56, KRK_TwoArea.lineInductance56, KRK_TwoArea.lineCapacitance56, KRK_TwoArea.lineConductance56);

  //Line56
	auto line56SP = SP::Ph1::PiLine::make("PiLine56", Logger::Level::debug);
	line56SP->setParameters(KRK_TwoArea.lineResistance56, KRK_TwoArea.lineInductance56, KRK_TwoArea.lineCapacitance56, KRK_TwoArea.lineConductance56);

  //Line67
	auto line67SP = SP::Ph1::PiLine::make("PiLine67", Logger::Level::debug);
	line67SP->setParameters(KRK_TwoArea.lineResistance67, KRK_TwoArea.lineInductance67, KRK_TwoArea.lineCapacitance67, KRK_TwoArea.lineConductance67);

  //Line78_1
	auto line78_1SP = SP::Ph1::PiLine::make("PiLine78_1", Logger::Level::debug);
	line78_1SP->setParameters(KRK_TwoArea.lineResistance78, KRK_TwoArea.lineInductance78, KRK_TwoArea.lineCapacitance78, KRK_TwoArea.lineConductance78);

  //Line78_2
	auto line78_2SP = SP::Ph1::PiLine::make("PiLine78_2", Logger::Level::debug);
	line78_2SP->setParameters(KRK_TwoArea.lineResistance78, KRK_TwoArea.lineInductance78, KRK_TwoArea.lineCapacitance78, KRK_TwoArea.lineConductance78);

  //Line89_1
	auto line89_1SP = SP::Ph1::PiLine::make("PiLine89_1", Logger::Level::debug);
	line89_1SP->setParameters(KRK_TwoArea.lineResistance89, KRK_TwoArea.lineInductance89, KRK_TwoArea.lineCapacitance89, KRK_TwoArea.lineConductance89);

  //Line89_2
	auto line89_2SP = SP::Ph1::PiLine::make("PiLine89_2", Logger::Level::debug);
	line89_2SP->setParameters(KRK_TwoArea.lineResistance89, KRK_TwoArea.lineInductance89, KRK_TwoArea.lineCapacitance89, KRK_TwoArea.lineConductance89);

  //Line910
	auto line910SP = SP::Ph1::PiLine::make("PiLine910", Logger::Level::debug);
	line910SP->setParameters(KRK_TwoArea.lineResistance910, KRK_TwoArea.lineInductance910, KRK_TwoArea.lineCapacitance910, KRK_TwoArea.lineConductance910);

  //Line1011
	auto line1011SP = SP::Ph1::PiLine::make("PiLine1011", Logger::Level::debug);
	line1011SP->setParameters(KRK_TwoArea.lineResistance1011, KRK_TwoArea.lineInductance1011, KRK_TwoArea.lineCapacitance1011, KRK_TwoArea.lineConductance1011);

	// Topology
	gen1SP->connect({ n1SP });
	gen2SP->connect({ n2SP });
	gen3SP->connect({ n3SP });
	gen4SP->connect({ n4SP });

	trafo15SP->connect({n1SP, n5SP});
	trafo26SP->connect({n2SP, n6SP});
	trafo311SP->connect({n3SP, n11SP});
	trafo410SP->connect({n4SP, n10SP});

	load7SP->connect({ n7SP });
	load9SP->connect({ n9SP });

	line56SP->connect({ n5SP, n6SP });
	line67SP->connect({ n6SP, n7SP });
	line78_1SP->connect({ n7SP, n8SP });
	line78_2SP->connect({ n7SP, n8SP });
	line89_1SP->connect({ n8SP, n9SP });
	line89_2SP->connect({ n8SP, n9SP });
	line910SP->connect({ n9SP, n10SP });
	line1011SP->connect({ n10SP, n11SP });

	auto systemSP = SystemTopology(60,
			SystemNodeList{n1SP, n2SP, n3SP, n4SP, n5SP, n6SP, n7SP, n8SP, n9SP, n10SP, n11SP},
			SystemComponentList{gen1SP, gen2SP, gen3SP, gen4SP, trafo15SP, trafo26SP, trafo311SP, trafo410SP, load7SP, load9SP, line56SP, line67SP, line78_1SP, line78_2SP, line89_1SP, line89_2SP, line910SP, line1011SP});

	// Initialization of dynamic topology
	systemSP.initWithPowerflow(systemPF, CPS::Domain::SP);

	// Decoupled Simulation
	IdentifiedObject::List components8_1;
	components8_1.push_back(line78_1SP);
	components8_1.push_back(line78_2SP);

	IdentifiedObject::List components8_2;
	components8_2.push_back(line89_1SP);
	components8_2.push_back(line89_2SP);

	decoupleNode(systemSP, "n8", components8_1, components8_2, delay, method);

	// Get pointers to added components
	auto n8_1SP = systemSP.node<SP::SimNode>("n8_1");
	auto n8_2SP = systemSP.node<SP::SimNode>("n8_2");
	line78_1SP = systemSP.component<SP::Ph1::PiLine>("PiLine78_1");
	line78_2SP = systemSP.component<SP::Ph1::PiLine>("PiLine78_2");
	line89_1SP = systemSP.component<SP::Ph1::PiLine>("PiLine89_1");
	line89_2SP = systemSP.component<SP::Ph1::PiLine>("PiLine89_2");

	// Logging
	auto loggerSP = DataLogger::make(simNameSP);
	loggerSP->logAttribute("v1", n1SP->attribute("v"));
	loggerSP->logAttribute("v2", n2SP->attribute("v"));
	loggerSP->logAttribute("v3", n3SP->attribute("v"));
	loggerSP->logAttribute("v4", n4SP->attribute("v"));
	loggerSP->logAttribute("v5", n5SP->attribute("v"));
	loggerSP->logAttribute("v6", n6SP->attribute("v"));
	loggerSP->logAttribute("v7", n7SP->attribute("v"));
	loggerSP->logAttribute("v8_1", n8_1SP->attribute("v"));
	loggerSP->logAttribute("v8_2", n8_2SP->attribute("v"));
	loggerSP->logAttribute("v9", n9SP->attribute("v"));
	loggerSP->logAttribute("v10", n10SP->attribute("v"));
	loggerSP->logAttribute("v11", n11SP->attribute("v"));
	loggerSP->logAttribute("v_line56", line56SP->attribute("v_intf"));
	loggerSP->logAttribute("i_line56", line56SP->attribute("i_intf"));
	loggerSP->logAttribute("v_line67", line67SP->attribute("v_intf"));
	loggerSP->logAttribute("i_line67", line67SP->attribute("i_intf"));
	loggerSP->logAttribute("v_line78_1", line78_1SP->attribute("v_intf"));
	loggerSP->logAttribute("i_line78_1", line78_1SP->attribute("i_intf"));
	loggerSP->logAttribute("v_line78_2", line78_2SP->attribute("v_intf"));
	loggerSP->logAttribute("i_line78_2", line78_2SP->attribute("i_intf"));
	loggerSP->logAttribute("v_line89_1", line89_1SP->attribute("v_intf"));
	loggerSP->logAttribute("i_line89_1", line89_1SP->attribute("i_intf"));
	loggerSP->logAttribute("v_line89_2", line89_2SP->attribute("v_intf"));
	loggerSP->logAttribute("i_line89_2", line89_2SP->attribute("i_intf"));
	loggerSP->logAttribute("v_line910", line910SP->attribute("v_intf"));
	loggerSP->logAttribute("i_line910", line910SP->attribute("i_intf"));
	loggerSP->logAttribute("v_line1011", line1011SP->attribute("v_intf"));
	loggerSP->logAttribute("i_line1011", line1011SP->attribute("i_intf"));
	loggerSP->logAttribute("v_gen1", gen1SP->attribute("v_intf"));
	loggerSP->logAttribute("i_gen1", gen1SP->attribute("i_intf"));
	loggerSP->logAttribute("v_gen2", gen2SP->attribute("v_intf"));
	loggerSP->logAttribute("i_gen2", gen2SP->attribute("i_intf"));
	loggerSP->logAttribute("v_gen3", gen3SP->attribute("v_intf"));
	loggerSP->logAttribute("i_gen3", gen3SP->attribute("i_intf"));
	loggerSP->logAttribute("v_gen4", gen4SP->attribute("v_intf"));
	loggerSP->logAttribute("i_gen4", gen4SP->attribute("i_intf"));
	loggerSP->logAttribute("v_load7", load7SP->attribute("v_intf"));
	loggerSP->logAttribute("i_load7", load7SP->attribute("i_intf"));
	loggerSP->logAttribute("v_load9", load9SP->attribute("v_intf"));
	loggerSP->logAttribute("i_load9", load9SP->attribute("i_intf"));

	// trafo
	loggerSP->logAttribute("i_trafo15", trafo15SP->attribute("i_intf"));
	loggerSP->logAttribute("v_trafo15", trafo15SP->attribute("v_intf"));
	loggerSP->logAttribute("i_trafo26", trafo26SP->attribute("i_intf"));
	loggerSP->logAttribute("v_trafo26", trafo26SP->attribute("v_intf"));
	loggerSP->logAttribute("i_trafo311", trafo311SP->attribute("i_intf"));
	loggerSP->logAttribute("v_trafo311", trafo311SP->attribute("v_intf"));
	loggerSP->logAttribute("i_trafo410", trafo410SP->attribute("i_intf"));
	loggerSP->logAttribute("v_trafo410", trafo410SP->attribute("v_intf"));

	auto itm = systemSP.component<Signal::DecouplingIdealTransformer_SP_Ph1>("itm_n8");
	loggerSP->logAttribute("v_itm", itm->attribute("v_ref"));
	loggerSP->logAttribute("i_itm", itm->attribute("i_ref"));

	Simulation simSP(simNameSP, Logger::Level::debug);
	simSP.setSystem(systemSP);
	simSP.setTimeStep(timeStep);
	simSP.setFinalTime(finalTime);
	simSP.setDomain(Domain::SP);
	simSP.doSplitSubnets(true);
  	simSP.doInitFromNodesAndTerminals(true);
	simSP.addLogger(loggerSP);

	simSP.run();
}

int main(int argc, char* argv[]) {


	//Simulation parameters
	String simName="SP_SynGenTrStab_KRK_TwoAreaTrafo_SteadyState_split_ITM";
	Real finalTime = 1.0;
	Real timeStep = 0.001;
	Bool startFaultEvent=true;
	Bool endFaultEvent=true;
	Real startTimeFault=10;
	Real endTimeFault=10.2;
	Real cmdInertia_G1 = KRK_TwoArea.H_G1;
	Real cmdInertia_G2 = KRK_TwoArea.H_G2;
  Real cmdInertia_G3 = KRK_TwoArea.H_G3;
  Real cmdInertia_G4 = KRK_TwoArea.H_G4;
	Real cmdDamping_G1=1.0;
	Real cmdDamping_G2=1.0;
  Real cmdDamping_G3=1.0;
  Real cmdDamping_G4=1.0;

	CommandLineArgs args(argc, argv);

  Real delay = 0.0001;
  String cosimMethodStr = "delay";
  String prefix = "1e-4";

	if (argc > 1) {
		timeStep = args.timeStep;
		finalTime = args.duration;
		if (args.name != "dpsim")
			simName = args.name;
		if (args.options.find("SCALEINERTIA_G1") != args.options.end())
			cmdInertia_G1 = args.getOptionReal("SCALEINERTIA_G1");
		if (args.options.find("SCALEINERTIA_G2") != args.options.end())
			cmdInertia_G2 = args.getOptionReal("SCALEINERTIA_G2");
    if (args.options.find("SCALEINERTIA_G3") != args.options.end())
			cmdInertia_G3 = args.getOptionReal("SCALEINERTIA_G3");
    if (args.options.find("SCALEINERTIA_G4") != args.options.end())
			cmdInertia_G4 = args.getOptionReal("SCALEINERTIA_G4");
		if (args.options.find("SCALEDAMPING_G1") != args.options.end())
			cmdDamping_G1 = args.getOptionReal("SCALEDAMPING_G1");
		if (args.options.find("SCALEDAMPING_G2") != args.options.end())
			cmdDamping_G2 = args.getOptionReal("SCALEDAMPING_G2");
    if (args.options.find("SCALEDAMPING_G3") != args.options.end())
			cmdDamping_G3 = args.getOptionReal("SCALEDAMPING_G3");
    if (args.options.find("SCALEDAMPING_G4") != args.options.end())
			cmdDamping_G4 = args.getOptionReal("SCALEDAMPING_G4");
		if (args.options.find("STARTTIMEFAULT") != args.options.end())
			startTimeFault = args.getOptionReal("STARTTIMEFAULT");
		if (args.options.find("ENDTIMEFAULT") != args.options.end())
			endTimeFault = args.getOptionReal("ENDTIMEFAULT");
    if (args.options.find("delay") != args.options.end())
      delay = args.getOptionReal("delay");
    if (args.options.find("method") != args.options.end())
      cosimMethodStr = args.getOptionString("method");
    if (args.options.find("prefix") != args.options.end())
      prefix = args.getOptionString("prefix");
	}

  String simNameDecoupled = simName + "_" + prefix;

  CouplingMethod cosimMethod = CouplingMethod::DELAY;
  if (cosimMethodStr == "extrapolation-zoh")
    cosimMethod = CouplingMethod::EXTRAPOLATION_ZOH;
  else if (cosimMethodStr == "extrapolation-linear")
    cosimMethod = CouplingMethod::EXTRAPOLATION_LINEAR;

	SP_SynGenTrStab_KRK_TwoAreaTrafo_SteadyState(simNameDecoupled, timeStep, finalTime, startFaultEvent, endFaultEvent, startTimeFault, endTimeFault, cmdInertia_G1, cmdInertia_G2, cmdInertia_G3, cmdInertia_G4, cmdDamping_G1, cmdDamping_G2, cmdDamping_G3, cmdDamping_G4, delay, cosimMethod);
}
