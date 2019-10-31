from sls_sim.SystemModel import LTISystem
from sls_sim.Simulator import Simulator
from sls_sim.SynthesisAlgorithm import SLS
from sls_sim.NoiseModel import *
from sls_sim.PlantGenerator import *
import numpy as np

def state_fdbk_example():
    sys = LTISystem (
        Nx = 10, Nw = 10
    )

    # generate sys._A, sys._B2
    GenerateDoublyStochasticChain (
        system_model = sys,
        rho = 1,
        actuator_density = 1,
        alpha = 0.2
    )

    # specify system matrices
    sys._B1  = np.eye (sys._Nx)
    sys._C1  = np.concatenate ((np.eye(sys._Nx), np.zeros([sys._Nu, sys._Nx])), axis = 0)
    sys._D12 = np.concatenate ((np.zeros([sys._Nx, sys._Nu]), np.eye(sys._Nu)), axis = 0)

    sim_horizon = 25
    simulator = Simulator (
        system = sys,
        horizon = sim_horizon
    )

    noise = FixedNoiseVector (Nw = sys._Nx, horizon = sim_horizon)
    noise.generateNoiseFromNoiseModel (cls = ZeroNoise)
    noise._w[0][sys._Nx/2] = 10

    sys.useNoiseModel (noise_model = noise)

    ## (1) basic sls (centralized controller)
    synthesizer = SLS (
        FIR_horizon = 20,
        obj_type = SLS.Objective.H2
    )

    synthesizer.setSystemModel (sys)

    controller = synthesizer.synthesizeControllerModel ()

    simulator.setController (controller=controller)

    sys.initialize (x0 = np.zeros([sys._Nx, 1]))
    controller.initialize ()
    x,y,z,u = simulator.run ()

if __name__ == '__main__':
    state_fdbk_example()