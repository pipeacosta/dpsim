import sys
import os.path
import logging
import json
import time

from datetime import datetime
from villas.node.node import Node as VILLASnode

from multiprocessing import Process

import dpsimpy
import dpsimpyvillas

base = os.path.splitext(os.path.basename(sys.argv[0]))[0]
log = logging.getLogger(base)

def dpsim0():
    sim_name = "FileImport0"
    time_step = 0.001
    final_time = 10

    n1 = dpsimpy.dp.SimNode('n1', dpsimpy.PhaseType.Single, [10])
    n2 = dpsimpy.dp.SimNode('n2', dpsimpy.PhaseType.Single, [5])

    evs = dpsimpy.dp.ph1.VoltageSource('v_intf', dpsimpy.LogLevel.debug)
    evs.set_parameters(complex(5, 0))

    vs1 = dpsimpy.dp.ph1.VoltageSource('vs_1', dpsimpy.LogLevel.debug)
    vs1.set_parameters(complex(10, 0))

    r12 = dpsimpy.dp.ph1.Resistor('r_12', dpsimpy.LogLevel.debug)
    r12.set_parameters(1)

    evs.connect([dpsimpy.dp.SimNode.gnd, n2])
    vs1.connect([dpsimpy.dp.SimNode.gnd, n1])
    r12.connect([n1, n2])

    sys = dpsimpy.SystemTopology(50, [n1, n2], [evs, vs1, r12])

    dpsimpy.Logger.set_log_dir('logs/' + sim_name)

    logger = dpsimpy.Logger(sim_name)
    logger.log_attribute('v1', 'v', n1)
    logger.log_attribute('v2', 'v', n2)
    logger.log_attribute('r12', 'i_intf', r12)
    logger.log_attribute('ievs', 'i_intf', evs)
    logger.log_attribute('vevs', 'v_intf', evs)

    sim = dpsimpy.RealTimeSimulation(sim_name)
    sim.set_system(sys)
    sim.set_time_step(time_step)
    sim.set_final_time(final_time)

    file_config = '''{
        "type": "file",
        "uri": "dpsim-mqtt-distributed-villas-results-22-09-12_16_27_53.csv",
        "in": {
            "epoch_mode": "original",
            "signals": [
                {
                    "name": "v_intf",
                    "type": "complex"
                }
            ]
        }
    }'''
    
    intf = dpsimpyvillas.InterfaceVillas(name="file-dpsim0", config=file_config)

    sim.add_interface(intf)
    sim.import_attribute(evs.attr('V_ref'), 0)

    return sim, intf

def test_file_import():
    logging.basicConfig(format='[%(asctime)s %(name)s %(levelname)s] %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)

    sim, intf = dpsim0() #intf needs to be extracted from the dpsim-function since the interface object gets deleted otherwise leading to SegFault when starting the simulation
    
    sim.run(1)

if __name__ == '__main__':
    test_file_import()