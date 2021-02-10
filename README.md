# System Level Synthesis
This repository includes the MATLAB (SLS-MATLAB) and Python (SLSpy) implementations of [System Level Synthesis](https://arxiv.org/abs/1904.01634). For questions, comments, and bug reports, please email the appropriate developer.

## SLS-MATLAB
SLS-MATLAB is maintained and developed by jsli (at) caltech (dot) edu. For more details, refer to the readme in the `sls-code/matlab` directory.

Bibtex entry for citation:

    @misc{Li2019_SLSMatlab,   
      title  = {{SLS-MATLAB}: Matlab Toolbox for System Level Synthesis},   
      author = {Li, Jing Shuang},   
      url    = {https://github.com/sls-caltech/sls-code},   
      year   = 2019   
    }

## SLSpy
SLSpy is maintained and developed by shtseng (at) caltech (dot) edu. For more details, refer to the README in the `sls-code/python` directory.

Bibtex entry for citation:

    @misc{Tseng2020_SLSpy,   
      title  = {{SLSpy}: Python-Based System-Level Controller Synthesis Framework},   
      author = {Tseng, Shih-Hao and Li, Jing Shuang},   
      url    = {http://arxiv.org/abs/2004.12565},   
      year   = {2020}   
    }   

SLSpy is included as a submodule. To access it, enter the 'sls-code/python' directory and run the following commands

    git submodule init
    git submodule update 

To enable automatic update of the submodule, run

    git config --global submodule.recurse true

This will enable automatic update whenever running a git pull.