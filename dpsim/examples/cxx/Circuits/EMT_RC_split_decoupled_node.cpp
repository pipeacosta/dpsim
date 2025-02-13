/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include "dpsim-models/Definitions.h"
#include "dpsim-models/EMT/EMT_Ph1_Capacitor.h"
#include "dpsim-models/EMT/EMT_Ph1_Resistor.h"
#include "dpsim-models/Signal/DecouplingIdealTransformer_EMT_Ph1.h"
#include "dpsim-models/SimNode.h"
#include <csignal>
#include <cstdlib>
#include <iostream>

#include <DPsim.h>
#include <dpsim/ThreadLevelScheduler.h>

using namespace DPsim;
using namespace CPS;

void decoupleNode(SystemTopology &sys, const String &nodeName, const IdentifiedObject::List &componentsAt1,
                  const IdentifiedObject::List &componentsAt2, Real ITMDelay, String method,
                  Matrix irLine_0, Real i_inf_0) {

  CouplingMethod cosimMethod = CouplingMethod::DELAY;

  if (method == "extrapolation-zoh")
    cosimMethod = CouplingMethod::EXTRAPOLATION_ZOH;
  else if (method == "extrapolation-linear")
    cosimMethod = CouplingMethod::EXTRAPOLATION_LINEAR;

  SimNode<Real>::List newNodes;
  SimPowerComp<Real>::List newComponents;

  auto intfNode = sys.node<EMT::SimNode>(nodeName);
  std::shared_ptr<TopologicalNode> nodeCopy1Topo = intfNode->clone(nodeName + "_1");
  std::shared_ptr<TopologicalNode> nodeCopy2Topo = intfNode->clone(nodeName + "_2");

  auto nodeCopy1 = std::dynamic_pointer_cast<SimNode<Real>>(nodeCopy1Topo);
  auto nodeCopy2 = std::dynamic_pointer_cast<SimNode<Real>>(nodeCopy2Topo);

  newNodes.push_back(nodeCopy1);
  newNodes.push_back(nodeCopy2);

  for (auto genComp : componentsAt1) {
    auto comp = std::dynamic_pointer_cast<SimPowerComp<Real>>(genComp);
    auto compCopy = comp->clone(comp->name());

    SimNode<Real>::List nodeCopies;
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
    auto comp = std::dynamic_pointer_cast<SimPowerComp<Real>>(genComp);
    auto compCopy = comp->clone(comp->name());

    SimNode<Real>::List nodeCopies;
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
  for (auto comp : newComponents)
    sys.addComponent(comp);

  auto idealTrafo = Signal::DecouplingIdealTransformer_EMT_Ph1::make("itm_" + nodeName,
                                                                Logger::Level::debug);
  idealTrafo->setParameters(nodeCopy1, nodeCopy2, ITMDelay, irLine_0, i_inf_0, cosimMethod);
  sys.addComponent(idealTrafo);
  sys.addComponents(idealTrafo->getComponents());
}

void doSim(String &name, SystemTopology &sys, Int threads, Real ts, bool isDecoupled = false) {

  // Logging
  auto logger = DataLogger::make(name);
  logger->logAttribute("v_0", sys.node<EMT::SimNode>("n0")->attribute("v"));
  logger->logAttribute("v_1", sys.node<EMT::SimNode>("n1")->attribute("v"));
	logger->logAttribute("i_rline", sys.component<EMT::Ph1::Resistor>("r_line")->mIntfCurrent, 1, 1);

  if (isDecoupled) {
    logger->logAttribute("v_2_1", sys.node<EMT::SimNode>("n2_1")->attribute("v"));
    logger->logAttribute("v_2_2", sys.node<EMT::SimNode>("n2_2")->attribute("v"));
    logger->logAttribute("i_intf", sys.component<Signal::DecouplingIdealTransformer_EMT_Ph1>("itm_n2")->attribute("i_intf"), 1, 1);
    logger->logAttribute("i_ref", sys.component<Signal::DecouplingIdealTransformer_EMT_Ph1>("itm_n2")->attribute("i_ref"), 1, 1);
    logger->logAttribute("v_intf", sys.component<Signal::DecouplingIdealTransformer_EMT_Ph1>("itm_n2")->attribute("v_intf"), 1, 1);
    logger->logAttribute("v_ref", sys.component<Signal::DecouplingIdealTransformer_EMT_Ph1>("itm_n2")->attribute("v_ref"), 1, 1);
  } else {
	  logger->logAttribute("v_2", sys.node<EMT::SimNode>("n2")->attribute("v"));
  }

  Simulation sim(name, Logger::Level::debug);
  sim.setSystem(sys);
  sim.setTimeStep(ts);
  sim.setFinalTime(1.2);
  sim.setDomain(Domain::EMT);
  sim.doSplitSubnets(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.addLogger(logger);
  if (threads > 0)
    sim.setScheduler(std::make_shared<OpenMPLevelScheduler>(threads));

  //std::ofstream of1("topology_graph.svg");
  //sys.topologyGraph().render(of1));

  sim.run();
  sim.logStepTimes(name + "_step_times");
}

SystemTopology buildTopology(String &name, float r1_r, float c1_c, float rLine_r, float r3_r,
                             float c2_c, Matrix n1_v0, Matrix n2_v0) {
  CPS::Logger::setLogDir("logs/" + name);

  // Nodes
	auto gnd = EMT::SimNode::GND;
  auto n0 = EMT::SimNode::make("n0");
	auto n1 = EMT::SimNode::make("n1");
	auto n2 = EMT::SimNode::make("n2");

	// Components
  auto vs = EMT::Ph1::VoltageSource::make("vs");
  vs->setParameters(Complex(1, 0), 50);
	auto r1 =  EMT::Ph1::Resistor::make("r_1");
	r1->setParameters(r1_r);
	auto c1 = EMT::Ph1::Capacitor::make("c_1");
	c1->setParameters(c1_c);
	auto rLine = EMT::Ph1::Resistor::make("r_line");
	rLine->setParameters(rLine_r);
	auto r3 =  EMT::Ph1::Resistor::make("r_3");
	r3->setParameters(r3_r);
	auto c2 = EMT::Ph1::Capacitor::make("c_2");
	c2->setParameters(c2_c);

	n1->setInitialVoltage(n1_v0 * PEAK1PH_TO_RMS3PH);
	n2->setInitialVoltage(n2_v0 * PEAK1PH_TO_RMS3PH);

	// Topology
  vs->connect({ gnd, n0 });
	r1->connect({ n0, n1 });
	rLine->connect({ n2, n1 });
	c1->connect({ n1, gnd });
	r3->connect({ n2, gnd });
	c2->connect({ n2, gnd });

	auto sys = SystemTopology(50,
		SystemNodeList{gnd, n0, n1, n2},
		SystemComponentList{vs,r1, c1, rLine, c2, r3});

	return sys;
}

int main(int argc, char *argv[]) {
  CommandLineArgs args(argc, argv);

  Int numThreads = 0;
  Int numSeq = 0;
  Real timeStep = 0.00005;
  Real delay = 0.0001;
  Real i_intf_0 = 0;
  String cosimMethod = "delay";
  String prefix = "1e-4";

  if (args.options.find("threads") != args.options.end())
    numThreads = args.getOptionInt("threads");
  if (args.timeStep)
    timeStep = args.timeStep;
  if (args.options.find("delay") != args.options.end())
    delay = args.getOptionReal("delay");
  if (args.options.find("i-intf-0") != args.options.end())
    i_intf_0 = args.getOptionReal("i-intf-0");
  if (args.options.find("method") != args.options.end())
    cosimMethod = args.getOptionString("method");
  if (args.options.find("prefix") != args.options.end())
    prefix = args.getOptionString("prefix");

  std::cout << "Simulate with " << numThreads
            << " threads, sequence number " << numSeq
            << ", co-simulation method " << cosimMethod << std::endl;

  float r1_r_1 = 0.1;
	float c1_c_1 = 1;
	float rLine_r_1 = 0.1;
	float r3_r_1 = 1;
	float c2_c_1 = 1;

  // Initial conditions, given by the problem
  Matrix n1_v0_1(1,1);
  n1_v0_1(0,0) = 0.0;
  Matrix n2_v0_1(1,1);
  n2_v0_1(0,0) = 0.0;

  Matrix irLine_0_1(1,1);
	irLine_0_1(0,0) = (n1_v0_1(0,0) - n2_v0_1(0,0)) / rLine_r_1;

  // Matrix ir3_0_1(1,1);
	// ir3_0_1(0,0) = (n2_v0_1(0,0)) / r3_r_1;

  // Monolithic Simulation
  String simNameMonolithic = "EMT_RC_monolithic";
  Logger::setLogDir("logs/" + simNameMonolithic);
  SystemTopology systemMonolithic = buildTopology(simNameMonolithic, r1_r_1, c1_c_1,
                                          rLine_r_1, r3_r_1, c2_c_1, n1_v0_1,
                                          n2_v0_1);

  doSim(simNameMonolithic, systemMonolithic, 0, timeStep);

  // Decoupled Simulation
  String simNameDecoupled = "EMT_RC_split_decoupled_" + prefix + "_" + std::to_string(numThreads) + "_" + std::to_string(numSeq);
  Logger::setLogDir("logs/" + simNameDecoupled);
  SystemTopology systemDecoupled = buildTopology(simNameDecoupled, r1_r_1, c1_c_1,
                                          rLine_r_1, r3_r_1, c2_c_1, n1_v0_1,
                                          n2_v0_1);

  IdentifiedObject::List components1;

  auto rLine = systemDecoupled.component<EMT::Ph1::Resistor>("r_line");
  components1.push_back(rLine);

  IdentifiedObject::List components2;
  auto c2 = systemDecoupled.component<EMT::Ph1::Capacitor>("c_2");
  components2.push_back(c2);
  // c2->setIntfVoltage(n2_v0_1);
  auto r3 = systemDecoupled.component<EMT::Ph1::Resistor>("r_3");
  // r3->setIntfCurrent(ir3_0_1);
  components2.push_back(r3);

  decoupleNode(systemDecoupled, "n2", components1, components2, delay, cosimMethod, irLine_0_1, i_intf_0);
  doSim(simNameDecoupled, systemDecoupled, numThreads, timeStep, true);
}
