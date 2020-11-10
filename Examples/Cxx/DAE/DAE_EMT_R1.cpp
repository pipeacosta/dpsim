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
using namespace CPS::EMT::Ph1;

int main(int argc, char* argv[])
{
	Real timeStep = 0.001;
	Real finalTime = 0.1;
	String simName = "DAE_EMT_VS_R1";
	Logger::setLogDir("logs/"+simName);

	// Nodes
	std::vector<Complex> initialVoltage_n1{ Complex(10,0), 0, 0 };
	auto n1 = SimNode::make("n1", PhaseType::Single, initialVoltage_n1);

	// Components
	auto vs = VoltageSource::make("vs", Logger::Level::info);
	vs->setParameters(Complex(10, 0), 50);
	auto r1 = Resistor::make("r_1");
	r1->setParameters(10);


	// Topology
	vs->connect(SimNode::List{ SimNode::GND, n1 });
	r1->connect(SimNode::List{ SimNode::GND, n1 });

	// Define system topology
	auto sys = SystemTopology(50, SystemNodeList{n1}, SystemComponentList{vs, r1});

	// Logger
	auto logger = DataLogger::make(simName);
	logger->addAttribute("vs", vs->attribute("v_intf"));
	logger->addAttribute("i_vs", vs->attribute("i_intf"));
	logger->addAttribute("vr", r1->attribute("v_intf"));
	logger->addAttribute("i_r", r1->attribute("i_intf"));
	logger->addAttribute("v1", n1->attribute("v"));

	Matrix vs_initCurrent = Matrix::Zero(1, 1);
	vs_initCurrent(0,0) = 1.0;
	vs->setIntfCurrent(vs_initCurrent);
	Simulation sim(simName, sys, timeStep, finalTime, Domain::EMT, Solver::Type::DAE);
	sim.doSplitSubnets(false);
	sim.addLogger(logger);
	sim.run();

	return 0;
}
