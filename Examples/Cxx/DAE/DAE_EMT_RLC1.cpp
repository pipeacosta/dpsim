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
	Real timeStep = 0.0001;
	Real finalTime = 0.1;
	String simName = "DAE_EMT_RLC1";
	Logger::setLogDir("logs/"+simName);

	// Nodes
	std::vector<Complex> initialVoltage{ 0.0, 0, 0 };
	auto n1 = SimNode::make("n1", PhaseType::Single, initialVoltage);
	auto n2 = SimNode::make("n2", PhaseType::Single, initialVoltage);
    auto n3 = SimNode::make("n3", PhaseType::Single, initialVoltage);

	// Components
	auto vs = VoltageSource::make("vs", Logger::Level::info);
	vs->setParameters(Complex(10, 0), 50);
	auto r1 = Resistor::make("r_1", Logger::Level::info);
	r1->setParameters(10);
	auto l1 = Inductor::make("l_1", Logger::Level::info);
	l1->setParameters(0.02);
    auto c1 = Capacitor::make("c_1", Logger::Level::info);
	c1->setParameters(0.0005);

	// Topology
	vs->connect(SimNode::List{ SimNode::GND, n1 });
	r1->connect(SimNode::List{ n1, n2 });
	l1->connect(SimNode::List{ n2, n3});
    c1->connect(SimNode::List{ n3, SimNode::GND});

	// Define system topology
	auto sys = SystemTopology(50, SystemNodeList{n1, n2, n3}, SystemComponentList{vs, r1, l1, c1});

	// Logger
	auto logger = DataLogger::make(simName);
	logger->addAttribute("i_c1", c1->attribute("i_intf"));
	logger->addAttribute("v_l", l1->attribute("v_intf"));
	logger->addAttribute("v_c", c1->attribute("v_intf"));
    logger->addAttribute("v_r", r1->attribute("v_intf"));

	Simulation sim(simName, sys, timeStep, finalTime, Domain::EMT, Solver::Type::DAE);
	sim.doSplitSubnets(false);
	sim.addLogger(logger);
	sim.run();

	return 0;
}
