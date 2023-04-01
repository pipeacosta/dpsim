import sys
import os.path
import logging

import dpsimpy
import villas.dataprocessing.timeseries as ts
import numpy as np

base = os.path.splitext(os.path.basename(sys.argv[0]))[0]
log = logging.getLogger(base)

def set_dpsim1(t_s, t_f, u_1_0):
    sim_name = "WeaklyCoupledCosim0"

    gnd = dpsimpy.emt.SimNode.gnd
    n1 = dpsimpy.emt.SimNode('n1', dpsimpy.PhaseType.Single, [5])
    n2 = dpsimpy.emt.SimNode('n2', dpsimpy.PhaseType.Single, [u_1_0])

    evs = dpsimpy.emt.ph1.VoltageSource('v_intf', dpsimpy.LogLevel.debug)
    evs.set_parameters(u_1_0)

    r_1 = dpsimpy.emt.ph1.Resistor("r_1")
    r_1.set_parameters(0.1)
    c_1 = dpsimpy.emt.ph1.Capacitor("c_1")
    c_1.set_parameters(1)
    r_line = dpsimpy.emt.ph1.Resistor('r_line', dpsimpy.LogLevel.debug)
    r_line.set_parameters(0.1)

    r_1.connect([n1, gnd])
    r_line.connect([n1, n2])
    c_1.connect([n1, gnd])
    evs.connect([gnd, n2])

    sys = dpsimpy.SystemTopology(50, [gnd, n1, n2], [evs, r_1, c_1, r_line])

    sim = dpsimpy.Simulation(sim_name)
    sim.set_domain(dpsimpy.Domain.EMT)
    sim.set_system(sys)
    sim.set_time_step(t_s)
    sim.set_final_time(t_f)
    # sim.do_steady_state_init(True)
    
    dpsimpy.Logger.set_log_dir('logs/' + sim_name)

    logger = dpsimpy.Logger(sim_name)
    logger.log_attribute('1_v_1', 'v', n1)
    logger.log_attribute('2_v_2', 'v', n2)
    logger.log_attribute('3_i_rline', 'i_intf', r_line)
    logger.log_attribute('4_i_evs', 'i_intf', evs)
    logger.log_attribute('5_v_evs', 'v_intf', evs)
    
    sim.add_logger(logger)
  
    return sim


def set_dpsim2(t_s, t_f, u_2_0):
    sim_name = "WeaklyCoupledCosim1"

    gnd = dpsimpy.emt.SimNode.gnd
    n2 = dpsimpy.emt.SimNode('n2', dpsimpy.PhaseType.Single, [2])

    ecs = dpsimpy.emt.ph1.CurrentSource('i_intf', dpsimpy.LogLevel.debug)
    ecs.set_parameters(u_2_0)
    c_2 = dpsimpy.emt.ph1.Capacitor("c_2")
    c_2.set_parameters(1)
    r_load = dpsimpy.emt.ph1.Resistor('r_load', dpsimpy.LogLevel.debug)
    r_load.set_parameters(1)

    ecs.connect([gnd, n2])
    c_2.connect([gnd, n2])
    r_load.connect([gnd, n2])

    sys = dpsimpy.SystemTopology(50, [gnd, n2], [ecs, c_2, r_load])

    sim = dpsimpy.Simulation(sim_name)
    sim.set_system(sys)
    sim.set_time_step(t_s)
    sim.set_final_time(t_f)

    dpsimpy.Logger.set_log_dir('logs/' + sim_name)

    logger = dpsimpy.Logger(sim_name)
    logger.log_attribute('1_v_2', 'v', n2)
    logger.log_attribute('2_i_rload', 'i_intf', r_load)
    logger.log_attribute('3_v_ecs', 'v_intf', ecs)
    logger.log_attribute('4_i_ecs', 'i_intf', ecs)

    sim.add_logger(logger)
  
    return sim

if __name__ == '__main__':
    logging.basicConfig(format='[%(asctime)s %(name)s %(levelname)s] %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)

    time_step = 1e-3
    final_time = 1.0
    t_k = 0.0
   
    # Communication y_2_0 -> S_1 and initialization of S_1
    y_2_0 = 2
    sim1 = set_dpsim1(time_step, final_time, y_2_0)
    sim1.start()
    y_1_0 = sim1.get_idobj_attr("v_intf", "i_intf").derive_coeff(0,0).get()
    print(y_1_0)
    k1_i = sim1.next()
    y_1_0 = sim1.get_idobj_attr("v_intf", "i_intf").derive_coeff(0,0).get()
    print(y_1_0)
    # y_1_0_ts = ts.TimeSeries("v2", t_k, y_1_0_dp)
    # y_1_0_emt = y_1_0_ts.frequency_shift(freq=50)

    # Communication y_1_0 -> S_2 and initialization of S_2
    # sim2 = set_dpsim2(time_step, final_time, y_1_0)
    # sim2.start()

    # Switch to S_2 and get initial output
    # y_2_0_dp = sim2.get_idobj_attr("n2", "v").derive_coeff(0,0).get()
    # y_2_0_ts = ts.TimeSeries("v2", t_k, y_2_0_dp)
    # y_2_0 = y_2_0_ts.frequency_shift(freq=50)

    # Communication y_1_0 -> S_2
    # u_1_0 = y_2_0
    # sim1.

    # Switch to S_2
    # y_1_0 = sim2.get_idobj_attr("n2", "v").derive_coeff(0,0).get()

    # k1_i = sim1.next()
    # k2_i = sim2.next()
    
    # while(k1_i <= final_time and k2_i <= final_time):
    #     k1_i = sim1.next()
    #     k2_i = sim2.next()
