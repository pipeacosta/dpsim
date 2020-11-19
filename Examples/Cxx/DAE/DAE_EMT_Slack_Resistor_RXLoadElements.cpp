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
using namespace CPS;

int main(int argc, char* argv[]) {
	// Parameters
	Real frequency = 50;
	Real Vnom = 110e3;
	Matrix rline_param = Matrix::Zero(3, 3);
	rline_param <<
		20., 0, 0,
		0, 20., 0,
		0, 0, 20.;
    Matrix rload_param = Matrix::Zero(3, 3);
	rload_param <<
		8066.66667, 0, 0,
		0, 8066.66667, 0,
		0, 0, 8066.66667;
    Matrix iload_param = Matrix::Zero(3, 3);
	iload_param <<
		25.677, 0, 0,
		0, 25.677, 0,
		0, 0, 25.677;

	// ----- DYNAMIC SIMULATION -----
	Real timeStepEMT  = 0.0001;
	Real finalTimeEMT = 0.1;
	String simNameEMT = "DAE_EMT_Slack_Resistor_RXLoadElements";
	Logger::setLogDir("logs/" + simNameEMT);

	// Nodes
	std::vector<Complex> initialVoltage_n1{ Complex(Vnom, 0), 
											Complex(Vnom, 0)*SHIFT_TO_PHASE_B,
											Complex(Vnom, 0)*SHIFT_TO_PHASE_C
										  };
	auto n1EMT = SimNode<Real>::make("n1EMT", PhaseType::ABC, initialVoltage_n1);
	std::vector<Complex> initialVoltage_n2{ Complex(89591.9,221.579)*PEAK1PH_TO_RMS3PH, 
											Complex(-44604.1,-77699.7)*PEAK1PH_TO_RMS3PH,
											Complex(-44987.9,77478.1)*PEAK1PH_TO_RMS3PH
										  };
	auto n2EMT = SimNode<Real>::make("n2EMT", PhaseType::ABC, initialVoltage_n2);

	// Slack
	auto slackEMT = EMT::Ph3::NetworkInjection::make("slackEMT", Logger::Level::debug);
	MatrixComp voltageRef_slackEMT = MatrixComp::Zero(3, 1);
	voltageRef_slackEMT(0, 0) = Complex(Vnom, 0);
	voltageRef_slackEMT(1, 0) = Complex(Vnom, 0)*SHIFT_TO_PHASE_B;
	voltageRef_slackEMT(2, 0) = Complex(Vnom, 0)*SHIFT_TO_PHASE_C;
	slackEMT->setParameters(voltageRef_slackEMT, frequency);

	// Resistor Line
	auto line_resistorEMT = EMT::Ph3::Resistor::make("line_resistorEMT", Logger::Level::debug);
	line_resistorEMT->setParameters(rline_param);

	// Resistor Load
	auto load_resistorEMT = EMT::Ph3::Resistor::make("load_resistorEMT", Logger::Level::info);
	load_resistorEMT->setParameters(rload_param);

	// inductor Load
	auto load_inductorEMT = EMT::Ph3::Inductor::make("load_inductorEMT", Logger::Level::debug);
	load_inductorEMT->setParameters(iload_param);

	// Topology
	slackEMT->connect({ n1EMT });
	line_resistorEMT->connect({ n2EMT, n1EMT });
	load_resistorEMT->connect({ SimNode<Real>::GND, n2EMT });
	load_inductorEMT->connect({ SimNode<Real>::GND, n2EMT });
	auto systemEMT = SystemTopology(frequency, SystemNodeList{n1EMT, n2EMT}, SystemComponentList{slackEMT, line_resistorEMT, load_resistorEMT, load_inductorEMT});

	// Logger
	auto loggerEMT = DataLogger::make(simNameEMT);
	loggerEMT->addAttribute("v1", n1EMT->attribute("v"));
	loggerEMT->addAttribute("v2", n2EMT->attribute("v"));
	loggerEMT->addAttribute("i_slack",  slackEMT->attribute("i_intf"));

    Matrix initCurrent_slack = Matrix::Zero(3, 1);
    initCurrent_slack(0, 0) = 11.133908;
	initCurrent_slack(1, 0) = -15.161624;
	initCurrent_slack(2, 0) = 4.027716;

	slackEMT->setIntfCurrent(initCurrent_slack);
	Simulation simEMT(simNameEMT, systemEMT, timeStepEMT, finalTimeEMT, Domain::EMT, Solver::Type::DAE, Logger::Level::debug);
	simEMT.doSplitSubnets(false);
	simEMT.addLogger(loggerEMT);
	simEMT.run();

	return 0;
}
