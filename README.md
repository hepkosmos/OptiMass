## OPTIMASS: a package for the minimization of kinematic mass functions with constraints
---
[![License: LGPL v2.1+](https://img.shields.io/badge/License-LGPL%20v2.1+-blue.svg)](https://www.gnu.org/licenses/lgpl.html)

OPTIMASS is a package for the minimization of kinematic mass functions with various kinematic constraints for 
general event topologies, providing users with an ability to obtain mass variables optimized for specific event topologies.








## Installation
---

#### 1) Requirements
  C++ compiler (checked for gcc = 5.4) 
  Python (python >= 2.7)
  ROOT with MINUIT2 activated [URL](https://root.cern.ch): Check with $root-config --has-minuit2
  Autotools and Libtools


#### 2) Downloading OPTIMASS

> $git clone https://github.com/hepkosmos/OptiMass.git


#### 3) Checking Relevant Operations by Interpreter

> $./optimass


#### 4) Building the Core Library of OptiMass

> $./optimass --build



## Generation of User's own Process Directory
---

### 1) Checking the List of the Process Cards (<myproc>.xml) in 'model_cards/' dir.

> $./optimass --list


### 6) Vim-editing a Process Card (<myproc>.xml) 

> $./optimass --vim <myproc>


### 7) Interpreting Users Process Cards and Generating a Job Directory:

> $./optimass --gen <myproc_1> <myproc_2> ... --dir <dir_path_name>


### 8) Working in the Job Directory:
Entered the process job directory, you can customize the `main.cpp` for loading your own events, 
simply it can be compiled by `make`, which generate the executable `optimass.x` 

'''
> $cd <dir_path_name> 
> $make
> $./optimass.x
'''


## Reference and Cite
---
When citing OptiMass, please use the following reference paper:

  OPTIMASS : A Package for the Minimization of Kinematic Mass Functions with Constraints
  [JHEP 1601(2016) 026](https://link.springer.com/article/10.1007%2FJHEP01%282016%29026) [arXiv:1508.00589 [hep-ph]](https://arxiv.org/abs/1508.00589v2)


