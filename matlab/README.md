# System Level Synthesis
MATLAB implementation of system level synthesis.

## Dependencies
This toolbox uses cvx, which can be downloaded from [here](http://cvxr.com/cvx/download/).

## Setup
Add the sls directory to MATLAB's directory by running `code/init.m`.

## Code
Reference materials

For `sls_base`
[System Level Synthesis, Anderson et al.](https://arxiv.org/pdf/1904.01634v1.pdf)
Equation numbers in the in-code documentation refer to equations from this paper unless otherwise specified.

For `sls_base/refine_controller.m`
[Separating Controller Design from Closed Loop Design, Li & Ho](https://arxiv.org/pdf/2006.05040.pdf)

For `sls_mpc`
[Distributed and Localized MPC via System Level Synthesis, Amo Alonso & Matni](https://arxiv.org/pdf/1909.10074.pdf)

### Key functions
`sls_base/state_fdbk_sls.m` contains the main SLS algorithm
`sls_base/refine_controller.m` implements the 2nd step of two-step SLS
`sls_base/simulate_state_fdbk.m` simulates the system using the SLS controller

`sls_mpc/sls_mpc.m` contains the main distributed MPC algorithm

### Examples
Suggested examples to start:
`sls_base/examples/state_fdbk_example.m`
`sls_mpc/examples/tests_algorithm_1.m`

### Classes
`LTISystem` contains all system matrices of the plant

`CLMaps` contains closed loop maps for the system

`Ctrller` contains controller matrices for the system (note that in conventional SLS, the CLMaps are directly used as controller matrices. To convert CLMaps directly to Ctrller, use Ctrller.ctrller_from_cl_maps(clMaps)).

`SimParams` contains parameters for simulation, such as disturbance

`SLSParams` contains parameters for the SLS algorithm

`MPCParams` contains parameters for the MPC algorithm