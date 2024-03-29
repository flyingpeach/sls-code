# SLS-MATLAB

## Dependencies & Setup
This toolbox uses cvx, which can be downloaded from [here](http://cvxr.com/cvx/download/). I strongly suggest that the user use Gurobi instead of SDPT3, as it is substantially faster.

The MATLAB code does not require installation. Simply add the sls directory and its subfolders to MATLAB's path; this can be done by opening MATLAB from the sls-code/matlab directory and running `init.m`.


## Standard SLS
Standard SLS is introduced in [System Level Synthesis](https://arxiv.org/abs/1904.01634).

`sls_base/state_fdbk_sls.m` contains the main algorithm for standard SLS
`sls_base/simulate_state_fdbk.m` simulates the system using this SLS controller

Suggested example file to try: `sls_base/examples/state_fdbk_example.m`


## Controller refinement
Controller refinement is a modification of standard SLS. It is run as a second step, after the standard SLS algorithm is run. It is introduced in [Separating Controller Design from Closed Loop Design](https://arxiv.org/abs/2006.05040). 

`sls_base/refine_controller.m` implements the controller refinement algorithm

Suggested example file to try: `sls_base/examples/ctrller_refinement_example.m`


## Blended SLS
Saturating states and/or inputs can be accommodated by blending multiple SLS controllers. This is described in [Achieving Performance and Safety in Large Scale Systems with Saturation using a Nonlinear System Level Synthesis Approach](https://arxiv.org/abs/2006.12766)

`sls_blend/two_blend_sls.m` contains an algorithm to blend two SLS controllers to accommodate saturation

`sls_blend/two_blend_sls.m` simulates the system using the blended controller

Suggested example file to try: `sls_blend/two_blend_example.m`


## Distributed MPC
Using the SLS parametrization, a distributed MPC algorithm can be developed. To the best of our knowledge, this is the first distributed MPC algorithm that is distributed in both synthesis and implementation, and is able to handle coupled objectives and constraints. The algorithm is introduced in [Distributed and Localized MPC via System Level Synthesis](https://arxiv.org/abs/1909.10074). This algorithm is made more computationally efficient in [Explicit Distributed and Localized Model Predictive Control via System Level Synthesis](https://arxiv.org/abs/2005.13807). The two-part paper Distributed and Localized Model Predictive Control [Part I: Synthesis and Implementation](https://arxiv.org/abs/2110.07010) and [Part II: Theoretical Guarantees](https://arxiv.org/abs/2203.00780) provides robust formulations and algorithms to calculate terminal set and cost.

`sls_mpc/mpc_uncoupled_distributed.m` and `sls_mpc/mpc_coupled_distributed.m` contain algorithms for nominal distributed MPC; `sls_mpc/rmpc_distributed.m` contains algorithms for robust distributed MPC.

Suggested example files to try: `sls_mpc/examples/nominal_uncoupled.m` and `sls_mpc/examples/nominal_coupled.m`


## D-Phi Iteration
Based on the SLS parametrization, the D-Phi iteration algorithm provides distributed, localized, and scalable robust control of systems with structured uncertainties. The algorithm is introduced in [Distributed Robust Control for Systems with Structured Uncertainties](https://arxiv.org/abs/2204.02493), and works with L1, L-infinity, and [nu](https://arxiv.org/abs/2204.05359) robustness.

`d_phi_iteration/d_phi_iteration.m` contains the algorithm that produces a distributed controller that is robust to structured uncertainty.

Suggested example files to try: `d_phi_iteration/d_phi_iteration_example.m`, which reproduces the simulation results of the paper
