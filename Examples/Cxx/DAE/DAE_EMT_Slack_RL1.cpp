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
    Real Vnom = 10;
    Real frequency = 50.0;
    Real inductance= 0.005;
    Real resistance = 5;
	Real timeStep = 0.001;
	Real finalTime = 0.1;
	String simName = "DAE_EMT_Slack_RL1";
	Logger::setLogDir("logs/" + simName);

	// Slack
	auto slack = NetworkInjection::make("slack", Logger::Level::info);
	MatrixComp voltageRef_slack = MatrixComp::Zero(3, 1);
	voltageRef_slack(0, 0) = Complex(Vnom, 0);
	voltageRef_slack(1, 0) = Complex(Vnom, 0)*SHIFT_TO_PHASE_B;
	voltageRef_slack(2, 0) = Complex(Vnom, 0)*SHIFT_TO_PHASE_C;
	slack->setParameters(voltageRef_slack, 50);

	// Inductor
	auto inductor1 = Inductor::make("Inductor", Logger::Level::info);
	Matrix ind_param = Matrix::Zero(3, 3);
	ind_param <<
		inductance, 0, 0,
		0, inductance, 0,
		0, 0, inductance;
	inductor1->setParameters(ind_param);

    // Resistor
	auto resistor1 = Resistor::make("Resistor", Logger::Level::info);
	Matrix resistor_param = Matrix::Zero(3, 3);
	resistor_param <<
		resistance, 0, 0,
		0, resistance, 0,
		0, 0, resistance;
	resistor1->setParameters(resistor_param);

    // Nodes
	std::vector<Complex> initialVoltage_n1{ voltageRef_slack(0, 0), 
											voltageRef_slack(1, 0),
											voltageRef_slack(2, 0)
										  };
	auto n1 = SimNode::make("n1", PhaseType::ABC, initialVoltage_n1);
	MatrixComp initCurrent_slack = MatrixComp::Zero(3, 1);
	Complex initCurrent = Complex(1.820339675, -0.571876575);
	Complex init_voltage_v2 = initCurrent*Complex(0, 2*PI*frequency*inductance);
    std::vector<Complex> initialVoltage_n2{ init_voltage_v2, 
											init_voltage_v2*SHIFT_TO_PHASE_B,
											init_voltage_v2*SHIFT_TO_PHASE_C
										  };
	auto n2 = SimNode::make("n2", PhaseType::ABC, initialVoltage_n2);

	// Topology
	slack->connect(SimNode::List{ n1 });
	resistor1->connect(SimNode::List{ n2, n1 });
    inductor1->connect(SimNode::List{ SimNode::GND, n2 });

	// Define system topology
	auto sys = SystemTopology(frequency, SystemNodeList{n1, n2}, SystemComponentList{slack, resistor1, inductor1});

	// Logger
	auto logger = DataLogger::make(simName);
	logger->addAttribute("v1", n1->attribute("v"));
    logger->addAttribute("v2", n2->attribute("v"));
    logger->addAttribute("v_slack", slack->attribute("v_intf"));
	logger->addAttribute("i_slack", slack->attribute("i_intf"));
    logger->addAttribute("v_resistor", resistor1->attribute("v_intf"));
	logger->addAttribute("i_resistor", resistor1->attribute("i_intf"));
	logger->addAttribute("v_inductor", inductor1->attribute("v_intf"));
	logger->addAttribute("i_inductor", inductor1->attribute("i_intf"));

	
	initCurrent_slack(0, 0) = initCurrent*RMS3PH_TO_PEAK1PH;
	initCurrent_slack(1, 0) = initCurrent_slack(0, 0)*SHIFT_TO_PHASE_B;
	initCurrent_slack(2, 0) = initCurrent_slack(0, 0)*SHIFT_TO_PHASE_C;

	slack->setIntfCurrent(initCurrent_slack.real());
	Simulation sim(simName, sys, timeStep, finalTime, Domain::EMT, Solver::Type::DAE);
	sim.doSplitSubnets(false);
	sim.addLogger(logger);
	sim.run();

	return 0;
}
