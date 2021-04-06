## OPTIMASS: a package for the minimization of kinematic mass functions with constraints
[![License: LGPL v2.1+](https://img.shields.io/badge/License-LGPL%20v2.1+-blue.svg)](https://www.gnu.org/licenses/lgpl.html)

OPTIMASS is a package for the minimization of kinematic mass functions with various kinematic constraints for 
general event topologies, providing users with an ability to obtain mass variables optimized for specific event topologies.








## Installation

#### 1) Requirements

1. C++ compiler (checked for gcc = 5.4) 
2. Python (python >= 2.7)
3. [ROOT](https://root.cern.ch) with [MINUIT2](https://seal.web.cern.ch/seal/MathLibs/Minuit2/html/) activated: check activation with `$root-config --has-minuit2`.
4. Autotools and Libtools


#### 2) Downloading OPTIMASS

    $git clone https://github.com/hepkosmos/OptiMass.git


#### 3) Checking relevant operations by interpreter

    $./optimass


#### 4) Building the core library of OPTIMASS

    $./optimass --build



## Generation of user's own process directory

#### 1) Checking the list of process cards `<myproc>.xml` in `'model_cards/'` directory:

    $./optimass --list


#### 2) Vim-editing a process card `<myproc>.xml`:

    $./optimass --vim <myproc>


#### 3) Generating a job directory for a set of processes:

    $./optimass --gen <myproc_1> <myproc_2> ... --dir <dir_path_name>


#### 4) Working in the job directory:
Entered the process job directory, you can customize the `main.cpp` for loading your own events, 
simply it can be compiled by `make`, which generate the executable `optimass.x` 

    $cd <dir_path_name> 
    $make
    $./optimass.x


## Reference and cite
When use OPTIMASS, please cite the following reference paper:

*OPTIMASS : A Package for the Minimization of Kinematic Mass Functions with Constraints 
[JHEP 1601(2016) 026](https://link.springer.com/article/10.1007%2FJHEP01%282016%29026) [arXiv:1508.00589 [hep-ph]](https://arxiv.org/abs/1508.00589v2)


