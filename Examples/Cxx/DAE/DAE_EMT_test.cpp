/**
 *
 * @author Markus Mirz <mmirz@eonerc.rwth-aachen.de>
 * @copyright 2017-2018, Institute for Automation of Complex Power Systems, EONERC
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
	String simName = "DAE_DP_test";
	Real timeStep = 0.00005;
	simName = "DAE_DP_test" + std::to_string(timeStep);
	Logger::setLogDir("logs/"+simName);

	// Nodes
	auto n1 = SimNode::make("n1");
	auto n2 = SimNode::make("n2");
    //auto n3 = Node::make("n3");

	// Components
	auto vs = VoltageSource::make("v_s");
	auto rl = Resistor::make("r_line");
    //auto ll = Inductor::make("l_line");
	auto rL = Resistor::make("r_load");

	// Topology
	vs->connect({ SimNode::GND, n1 });
	rl->connect({ n1, n2 });
    //ll->connect({ n2, n3 });
    rL->connect({ SimNode::GND, n2 });

	// Parameters
	vs->setParameters(Complex(10000, 0));
	rl->setParameters(1);
    //ll->setParameters(1);
	rL->setParameters(1000);

	auto sys = SystemTopology(50, SystemNodeList{SimNode::GND, n1, n2}, SystemComponentList{vs, rl, rL});
	Simulation sim(simName, sys, timeStep, 1.0, Domain::EMT, Solver::Type::DAE);
	sim.doSplitSubnets(false);
	sim.run();

	return 0;
}
