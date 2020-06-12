# System Level Synthesis
This repository includes the MATLAB (SLS-MATLAB) and Python (SLSpy) implementations of System Level Synthesis.

## How to access SLSpy
SLSpy is included as a submodule. To access it, enter the 'sls-code/python' directory and run the following commands

    git submodule init
    git submodule update 

To enable automatic update of the submodule, run

    git config --global submodule.recurse true

This will enable automatic update whenever running a git pull.