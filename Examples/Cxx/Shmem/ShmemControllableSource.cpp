/** Example of shared memory interface
 *
 * @author Markus Mirz <mmirz@eonerc.rwth-aachen.de>
 * @copyright 2017, Institute for Automation of Complex Power Systems, EONERC
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

#include "DPsim.h"

using namespace DPsim;
using namespace CPS::Components::DP;

int main(int argc, char *argv[]) {

	struct shmem_conf conf;
	conf.samplelen = 64;
	conf.queuelen = 1024;
	conf.polling = false;

	String in  = "/dpsim10";
	String out = "/dpsim01";

	ShmemInterface shmem(in, out, &conf);

	Real timeStep = 0.000150;

	// Nodes
	auto n1 = Node::make("n1");

	// Components
	auto ecs = CurrentSource::make("v_intf", Node::List{GND, n1}, Complex(100, 0));
	auto r1 = Resistor::make("r_1", Node::List{GND, n1}, 1);

	shmem.registerControlledAttribute(ecs->findAttribute<Complex>("current_ref"), 0, 1);
	shmem.registerExportedAttribute(ecs->findAttribute<Complex>("comp_current"), 0, 1);

	auto sys = SystemTopology(50, Node::List{GND, n1}, ComponentBase::List{ecs, r1});
	auto sim = Simulation("ShmemControllableSource", sys, timeStep, 1);

	sim.addInterface(&shmem);
	sim.run(false);

	return 0;
}