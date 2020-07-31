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

int main(int argc, char* argv[])
{
	Real timeStep = 0.0002;
	Real finalTime = 1;
	String simName = "DAE_EMT_Slack_PiLine_RXLoad";
	Logger::setLogDir("logs/" + simName);

	// Slack
	auto slack = NetworkInjection::make("slack", Logger::Level::info);
	MatrixComp initialVoltage_slack = MatrixComp::Zero(3, 1);
	initialVoltage_slack(0, 0) = Complex(110000, 0);
	initialVoltage_slack(1, 0) = initialVoltage_slack(0, 0) * SHIFT_TO_PHASE_B;
	initialVoltage_slack(2, 0) = initialVoltage_slack(0, 0) * SHIFT_TO_PHASE_C;
	slack->setParameters(initialVoltage_slack, 50);

	// PiLine
	auto PiLine = PiLine::make("PiLine", Logger::Level::info);
	Matrix rline_param = Matrix::Zero(3, 3);
	rline_param <<
		4.14, 0, 0,
		0, 4.14, 0,
		0, 0, 4.14;
	Matrix lf_param = Matrix(3, 3);
	lf_param <<
		0.0077031, 0, 0,
		0, 0.0077031, 0,
		0, 0, 0.0077031;
	Matrix cf_param = Matrix(3, 3);
	cf_param <<
		2.54648e-6, 0, 0,
		0, 2.54648e-6, 0,
		0, 0, 2.54648e-6;
	PiLine->setParameters(rline_param, lf_param, cf_param);

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
	std::vector<Complex> initialVoltage_n1{ initialVoltage_slack(0, 0), 
											initialVoltage_slack(1, 0), 
											initialVoltage_slack(2, 0)
										  };
	std::vector<Complex> initialVoltage_n2{ Complex(110000, 0), Complex(110000, 0)*SHIFT_TO_PHASE_B, Complex(110000, 0)*SHIFT_TO_PHASE_C };
	auto n1 = SimNode::make("n1", PhaseType::ABC, initialVoltage_n1);
	auto n2 = SimNode::make("n2", PhaseType::ABC, initialVoltage_n2);

	// Topology
	slack->connect(SimNode::List{ n1 });
	PiLine->connect(SimNode::List{ n1, n2 });
	rxLoad->connect(SimNode::List{ n2 });

	// Define system topology
	auto sys = SystemTopology(50, SystemNodeList{n1, n2}, SystemComponentList{slack, PiLine, rxLoad});

	// Logger
	auto logger = DataLogger::make(simName);
    // logger->addAttribute("vs",  vs->attribute("v_intf"));
	// logger->addAttribute("i_c", c1->attribute("i_intf"));
	// logger->addAttribute("v_l", l1->attribute("v_intf"));
	// logger->addAttribute("v_c", c1->attribute("v_intf"));
    // logger->addAttribute("v_r", r1->attribute("v_intf"));

	MatrixComp initCurrent_slack = MatrixComp::Zero(3, 1);
	initCurrent_slack(0, 0) = Complex(0, 0);
	initCurrent_slack(1, 0) = Complex(0, 0);
	initCurrent_slack(2, 0) = Complex(0, 0);
	slack->setInitialCurrent(initCurrent_slack);
	Simulation sim(simName, sys, timeStep, finalTime, Domain::EMT, Solver::Type::DAE);
	sim.doSplitSubnets(false);
	sim.addLogger(logger);
	sim.run();

	return 0;
}
