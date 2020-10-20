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

#define SHIFT_TO_PHASE_B Complex(cos(-2 * M_PI / 3), sin(-2 * M_PI / 3))
#define SHIFT_TO_PHASE_C Complex(cos(2 * M_PI / 3), sin(2 * M_PI / 3))
#define RMS3PH_TO_PEAK1PH sqrt(2./3.)
#define PEAK1PH_TO_RMS3PH sqrt(3./2.)
#define SQRT2 sqrt(2.)

int main(int argc, char* argv[])
{
	Real timeStep = 0.0001;
	Real finalTime = 0.1;
	String simName = "DAE_EMT_Slack_Resistor_RXLoad";
	Logger::setLogDir("logs/" + simName);

	// Slack
	auto slack = NetworkInjection::make("slack", Logger::Level::info);
	MatrixComp voltageRef_slack = MatrixComp::Zero(3, 1);
	voltageRef_slack(0, 0) = Complex(110e3, 0);
	voltageRef_slack(1, 0) = Complex(110e3, 0)*SHIFT_TO_PHASE_B;
	voltageRef_slack(2, 0) = Complex(110e3, 0)*SHIFT_TO_PHASE_C;
	slack->setParameters(voltageRef_slack, 50);

	// Resistor
	auto resistor = Resistor::make("Resistor", Logger::Level::info);
	Matrix resistor_param = Matrix::Zero(3, 3);
	resistor_param <<
		20, 0, 0,
		0, 20, 0,
		0, 0, 20;
	resistor->setParameters(resistor_param);

	// RXLoad
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

	// Node1
	std::vector<Complex> initialVoltage_n1{ voltageRef_slack(0, 0), 
											voltageRef_slack(1, 0),
											voltageRef_slack(2, 0)
										  };
	auto n1 = SimNode::make("n1", PhaseType::ABC, initialVoltage_n1);
	
	// Node2
	Complex init_v2 = Complex(109863.6368, 136.0259582);
	std::vector<Complex> initialVoltage_n2{ init_v2, 
											init_v2*SHIFT_TO_PHASE_B,
											init_v2*SHIFT_TO_PHASE_C
										  };
	auto n2 = SimNode::make("n2", PhaseType::ABC, initialVoltage_n2);

	// Topology
	slack->connect(SimNode::List{ n1 });
	resistor->connect(SimNode::List{ n2, n1 });
	rxLoad->connect(SimNode::List{ n2 });

	// Define system topology
	auto sys = SystemTopology(50, SystemNodeList{n1, n2}, SystemComponentList{slack, resistor, rxLoad});

	// Logger
	auto logger = DataLogger::make(simName);
	logger->addAttribute("v1", n1->attribute("v"));
	logger->addAttribute("v2", n1->attribute("v"));
	logger->addAttribute("i_slack",  slack->attribute("i_intf"));

	MatrixComp initCurrent_slack = MatrixComp::Zero(3, 1);
	Complex initCurrent = Complex(5.567004844, -5.55323649);
	initCurrent_slack(0, 0) = initCurrent;
	initCurrent_slack(1, 0) = initCurrent*SHIFT_TO_PHASE_B;
	initCurrent_slack(2, 0) = initCurrent*SHIFT_TO_PHASE_C;

	slack->setInitialCurrent(initCurrent_slack);
	Simulation sim(simName, sys, timeStep, finalTime, Domain::EMT, Solver::Type::DAE);
	sim.doSplitSubnets(false);
	sim.addLogger(logger);
	sim.run();

	return 0;
}
