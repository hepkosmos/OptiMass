## OPTIMASS: a package for the minimization of kinematic mass functions with constraints
[![License: LGPL v2.1+](https://img.shields.io/badge/License-LGPL%20v2.1+-blue.svg)](https://www.gnu.org/licenses/lgpl.html)

OPTIMASS is a package for the minimization of kinematic mass functions with various kinematic constraints for 
general event topologies, providing users with an ability to obtain mass variables optimized for specific decay processes.

Starting from user's process card in `XML` format
which defining the process information as follows:

   1. Decay chains with particle elements at nodes and leaves,
   2. Mass function of the masses of the particles,
   3. Invisible particles to be optimized for the mass function,
   4. Kinematic constraint functions, 
   5. Assignment of independent PT-conserving chains 
   (for multiple decay chains from independent events),

the interpreter (Python) of OPTIMASS can generate a process job directory including 

1. process dictionary codes (C++)
2. main function code (C++) for loading events and running OPTIMASS

for a selective set of user processes all together.

OPTIMASS's constrained minimization of the mass function with respect to the invisible momenta, subject to the kinematic constrains, is implemented by the Augmented Lagrange Method, utilizing the libraries of ROOT with MINUIT2 for a series of unconstrained minimizations required.



## Installation 

#### 1) Requirements

1. C++ compiler (tested with g++ v5.4-7.4) 
2. Python 2.7
3. [ROOT](https://root.cern.ch) with [MINUIT2](https://seal.web.cern.ch/seal/MathLibs/Minuit2/html/) activated: check activation with `$root-config --has-minuit2`.
4. Autotools and Libtools


#### 2) Downloading OPTIMASS (ver 2.0)
```bash
    $git clone https://github.com/hepkosmos/OptiMass.git
```

#### 3) Checking relevant operations by interpreter
```bash
    $./optimass
```

#### 4) Building the core library of OPTIMASS
```bash
    $./optimass --build
```


## Generation of user's own process directory

#### 1) Checking the list of process cards `<myproc>.xml` in `'model_cards/'` directory:
```bash
    $./optimass --list
```

#### 2) Vim-editing a process card `<myproc>.xml`:
```bash
    $./optimass --vim <myproc>
```

#### 3) Generating a job directory for a set of processes:
```bash
    $./optimass --gen <myproc_1> <myproc_2> ... --dir <dir_path_name>
```

#### 4) Working in the job directory:
Entered the process job directory, you can customize the `main.cpp` for loading your own events, 
simply it can be compiled by `make` to generate the executable `optimass.x` 
```bash
    $cd <dir_path_name> 
    $make
    $./optimass.x
```

## Citation

When you cite OPTIMASS, please use the script below:

``` bibtex
@article{Cho_2016,
   title={OPTIMASS: a package for the minimization of kinematic mass functions with constraints},
   volume={2016},
   ISSN={1029-8479},
   url={http://dx.doi.org/10.1007/JHEP01(2016)026},
   DOI={10.1007/jhep01(2016)026},
   number={1},
   journal={Journal of High Energy Physics},
   publisher={Springer Science and Business Media LLC},
   author={Cho, Won Sang and Gainer, James S. and Kim, Doojin and Lim, Sung Hak and Matchev, Konstantin T. and Moortgat, Filip and Pape, Luc and Park, Myeonghun},
   year={2016},
   month={Jan}
}
```

## References

* W.S. Cho et al, OPTIMASS : A Package for the Minimization of Kinematic Mass Functions with Constraints, [JHEP 1601(2016) 026](https://link.springer.com/article/10.1007%2FJHEP01%282016%29026), [arXiv:1508.00589 [hep-ph]](https://arxiv.org/abs/1508.00589v2)
* W.S. Cho et al, On-shell constrained M2​ variables with applications to mass measurements and topology disambiguation, [JHEP 08 (2014) 070](https://doi.org/10.1007/JHEP08(2014)070), [arXiv:1401.1449](https://arxiv.org/abs/1401.1449v2).
* C.B. Park, YAM2: Yet another library for the M2 variables using sequential quadratic programming, [Comput. Phys. Commun. 264 (2021) 107967](https://doi.org/10.1016/j.cpc.2021.107967), [arXiv:2007.15537](https://arxiv.org/abs/2007.15537).
* J. Nocedal and S. Wright, [Numerical Optimization](https://link.springer.com/book/10.1007/978-0-387-40065-5), Springer, 2006.
