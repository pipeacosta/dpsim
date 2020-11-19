/**
 *
 * @author Martin Moraga <mmoraga@eonerc.rwth-aachen.de>
 * @copyright 2017-2020, Institute for Automation of Complex Power Systems, EONERC
 *
 * DPsim
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

#include <DPsim.h>

using namespace DPsim;
using namespace CPS::EMT;
using namespace CPS::EMT::Ph3;

int main(int argc, char* argv[])
{
	Real timeStep = 0.001;
	Real finalTime = 0.1;
	String simName = "DAE_EMT_Slack_RXLoad";
	Logger::setLogDir("logs/" + simName);

	// Slack
	auto slack = NetworkInjection::make("slack", Logger::Level::info);
	MatrixComp voltageRef_slack = MatrixComp::Zero(3, 1);
	voltageRef_slack(0, 0) = Complex(110e3, 0);
	voltageRef_slack(1, 0) = Complex(110e3, 0)*SHIFT_TO_PHASE_B;
	voltageRef_slack(2, 0) = Complex(110e3, 0)*SHIFT_TO_PHASE_C;
	slack->setParameters(voltageRef_slack, 50);
    MatrixComp initCurrent_slack = MatrixComp::Zero(3, 1);
    initCurrent_slack(0, 0) = Real(11.13404383);
	initCurrent_slack(1, 0) = Real(-15.209387);
	initCurrent_slack(2, 0) = Real(4.075343);

	// ZLoad
	Matrix p_load = Matrix::Zero(3, 3);
	p_load <<
		0.5e6, 0, 0,
		0, 0.5e6, 0,
		0, 0, 0.5e6;
	Matrix q_load = Matrix::Zero(3, 3);
	q_load <<
		0.5e6, 0, 0,
		0, 0.5e6, 0,
		0, 0, 0.5e6;
	auto rxLoad = RXLoad::make("rxLoad", p_load, q_load, 110e3, Logger::Level::info);

	// Nodes
    std::vector<Complex> initialVoltage_n1{ voltageRef_slack(0, 0), 
											voltageRef_slack(1, 0),
											voltageRef_slack(2, 0)
										  };
	auto n1 = SimNode::make("n1", PhaseType::ABC, initialVoltage_n1);

	// Topology
	slack->connect(SimNode::List{ n1 });
	rxLoad->connect(SimNode::List{ n1 });

	// Define system topology
	auto sys = SystemTopology(50, SystemNodeList{n1}, SystemComponentList{slack, rxLoad});

	// Logger
	auto logger = DataLogger::make(simName);
	// logger->addAttribute("v1", n1->attribute("v"));
    logger->addAttribute("v_slack",  slack->attribute("v_intf"));
	logger->addAttribute("i_slack",  slack->attribute("i_intf"));
	logger->addAttribute("v_load",  rxLoad->attribute("v_intf"));
	logger->addAttribute("i_load",  rxLoad->attribute("i_intf"));

	slack->setIntfCurrent(initCurrent_slack.real());
	Simulation sim(simName, sys, timeStep, finalTime, Domain::EMT, Solver::Type::DAE);
	sim.doSplitSubnets(false);
	sim.addLogger(logger);
	sim.run();

	return 0;
}
