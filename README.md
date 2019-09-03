# System Level Synthesis
MATLAB implementation of system level synthesis and regularization for design (RFD).

## Setup
Add the sls directory to MATLAB's directory either manually or by running `code/init.m`.

## Code
All equations referenced in the documentation are from [this tutorial paper](https://arxiv.org/pdf/1904.01634v1.pdf).

### Examples
Several example files can be found in `code/examples`. Some of these (especially `cdc_examples.m`) generate *many* graphs, so it's suggested to run them one section at a time. Suggestion: start with `visually_pleasing_example.m` and `state_fdbk_example.m`.

### Classes
`LTISystem` contains all system matrices of the plant

`SimParams` contains parameters for simulation (i.e. time, disturbance)

`SLSOutputs` contains outputs of the SLS solver

`SLSParams` contains parameters for the SLS solver

### Key functions
`code/simulate_system.m` simulates the system according to the controller synthesized by the SLS solver

`code/state_fdbk_sls.m` is the main SLS solver

## Old version
The pre-github version of the code can be found on the [SLS Wiki](
http://slswiki.cms.caltech.edu/index.php/SLS_Toolbox). It is also tagged as a release on the repository; use `git checkout v1.0` to access it.