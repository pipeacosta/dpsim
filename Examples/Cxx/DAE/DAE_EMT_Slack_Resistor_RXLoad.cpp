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
using namespace CPS::CIM;

int main(int argc, char* argv[]) {
	// Parameters
	Real frequency = 50;
	Real Vnom = 110e3;
	//Real Vnom = 110;
	Matrix rline_param = Matrix::Zero(3, 3);
	///*
	rline_param <<
		20., 0, 0,
		0, 20., 0,
		0, 0, 20.;
	//	*/
	/*
	rline_param <<
		2., 0, 0,
		0, 2., 0,
		0, 0, 2.;	
	*/
	Real pLoadNom = 0.5e6;
	Real qLoadNom = 0.5e6;
	//Real pLoadNom = 10000;
	//Real qLoadNom = 10000;
	Real r_load = std::pow(Vnom/sqrt(3), 2) * (1/pLoadNom);
	Real i_load = std::pow(Vnom/sqrt(3), 2) * (1/qLoadNom) / (2 * PI * frequency);
	Matrix p_load = Matrix::Zero(3, 3);
	p_load <<
		pLoadNom, 0, 0,
		0, pLoadNom, 0,
		0, 0, pLoadNom;
	Matrix q_load = Matrix::Zero(3, 3);
	q_load <<
		qLoadNom, 0, 0,
		0, qLoadNom, 0,
		0, 0, qLoadNom;
	Matrix rloadPF_param = Matrix::Zero(3, 3);
	rloadPF_param <<
		r_load, 0, 0,
		0, r_load, 0,
		0, 0, r_load;
	Matrix iload_param = Matrix::Zero(3, 3);
	iload_param <<
		i_load, 0, 0,
		0, i_load, 0,
		0, 0, i_load;


	// ----- POWERFLOW FOR INITIALIZATION -----
	Real timeStepPF = 0.001;
	Real finalTimePF = 0.1;
	String simNamePF = "SP_Ph3_Slack_RLine_RXLoad_Init";
	Logger::setLogDir("logs/" + simNamePF);

	// Nodes
	auto n1PF = SimNode<Complex>::make("n1PF", PhaseType::ABC);
	auto n2PF = SimNode<Complex>::make("n2PF", PhaseType::ABC);

	// Components
	// voltage source
	auto vsPF = SP::Ph3::VoltageSource::make("vsPF");
	vsPF->setParameters(Vnom*RMS3PH_TO_PEAK1PH);
    //vsPF->setParameters(Vnom);

	// RLine
	auto rlinePF = SP::Ph3::Resistor::make("rlinePF");
	rlinePF->setParameters(rline_param);

	// RXLoad
    auto rloadPF = SP::Ph3::Resistor::make("rloadPF");
	rloadPF->setParameters(rloadPF_param);
	auto iloadPF = SP::Ph3::Inductor::make("iloadPF");
	iloadPF->setParameters(iload_param);

	// Topology
	vsPF->connect({SimNode<Complex>::GND, n1PF});
	rlinePF->connect({n2PF, n1PF});
	rloadPF->connect({ SimNode<Complex>::GND, n2PF });
	iloadPF->connect({ SimNode<Complex>::GND, n2PF });

	auto systemPF  = SystemTopology(frequency, SystemNodeList{ n1PF, n2PF}, SystemComponentList{ vsPF, rlinePF, rloadPF, iloadPF});

	// Logging
	auto loggerPF  = DataLogger::make(simNamePF);
	loggerPF->addAttribute("v1", n1PF->attribute("v"));
	loggerPF ->addAttribute("v2", n2PF->attribute("v"));
	loggerPF ->addAttribute("i_slack", rlinePF->attribute("i_intf"));
	loggerPF ->addAttribute("i_line", rlinePF->attribute("i_intf"));

	Simulation simPF(simNamePF, Logger::Level::info);
	simPF.setSystem(systemPF);
	simPF.addLogger(loggerPF);
	simPF.setDomain(Domain::SP);
	simPF.setTimeStep(timeStepPF);
	simPF.setFinalTime(finalTimePF);
	simPF.run();

	// ----- DYNAMIC SIMULATION -----
	Real timeStepEMT  = 0.0001;
	Real finalTimeEMT = 0.1;
	String simNameEMT = "DAE_EMT_Slack_RLine_RXLoad_PF_Init";
	Logger::setLogDir("logs/" + simNameEMT);

	// Nodes
	std::vector<Complex> initialVoltage_n1{ n1PF->voltage()(0,0).real()*PEAK1PH_TO_RMS3PH, 
											n1PF->voltage()(1,0).real()*PEAK1PH_TO_RMS3PH,
											n1PF->voltage()(2,0).real()*PEAK1PH_TO_RMS3PH
										  };
	auto n1EMT = SimNode<Real>::make("n1EMT", PhaseType::ABC, initialVoltage_n1);
	std::vector<Complex> initialVoltage_n2{ n2PF->voltage()(0,0).real()*PEAK1PH_TO_RMS3PH, 
											n2PF->voltage()(1,0).real()*PEAK1PH_TO_RMS3PH,
											n2PF->voltage()(2,0).real()*PEAK1PH_TO_RMS3PH
										  };
	auto n2EMT = SimNode<Real>::make("n2EMT", PhaseType::ABC, initialVoltage_n2);

	// Slack
	auto slackEMT = EMT::Ph3::NetworkInjection::make("slackEMT", Logger::Level::info);
	MatrixComp voltageRef_slackEMT = MatrixComp::Zero(3, 1);
	voltageRef_slackEMT(0, 0) = Complex(Vnom, 0);
	voltageRef_slackEMT(1, 0) = Complex(Vnom, 0)*SHIFT_TO_PHASE_B;
	voltageRef_slackEMT(2, 0) = Complex(Vnom, 0)*SHIFT_TO_PHASE_C;
	slackEMT->setParameters(voltageRef_slackEMT, frequency);

	// Resistor
	auto resistorEMT = EMT::Ph3::Resistor::make("ResistorEMT", Logger::Level::info);
	resistorEMT->setParameters(rline_param);

	// RXLoad
	auto rxLoadEMT = EMT::Ph3::RXLoad::make("rxLoadEMT", p_load, q_load, Vnom, Logger::Level::info);

	// Topology
	slackEMT->connect({ n1EMT });
	resistorEMT->connect({ n2EMT, n1EMT });
	rxLoadEMT->connect({ n2EMT });
	auto systemEMT = SystemTopology(frequency, SystemNodeList{n1EMT, n2EMT}, SystemComponentList{slackEMT, resistorEMT, rxLoadEMT});

	// Logger
	auto loggerEMT = DataLogger::make(simNameEMT);
	loggerEMT->addAttribute("v1", n1EMT->attribute("v"));
	loggerEMT->addAttribute("v2", n2EMT->attribute("v"));
	loggerEMT->addAttribute("i_slack",  slackEMT->attribute("i_intf"));

	slackEMT->setIntfCurrent(rlinePF->intfCurrent().real());
	Simulation simEMT(simNameEMT, systemEMT, timeStepEMT, finalTimeEMT, Domain::EMT, Solver::Type::DAE);
	simEMT.doSplitSubnets(false);
	simEMT.addLogger(loggerEMT);
	simEMT.run();

	return 0;
}
