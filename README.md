# ctpToMmoran

This C++ code was used to produce the Moran simulations discussed in the paper T. Lenaerts, J.M. Pacheco and F.C. Santos (2022) Evolution of a theory of mind. 
The parameters can be set in the main file, needing thus recompilation for each result.

## Contents

- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Compilation Guide](#compilation-guide)
- [Demo Guide](#demo)
- [Results](#demo)
- [Citation](#citation)
# Repo Contents

All *.CPP and *.HPP files:
- `rangen.h/rangen.cpp`
- `ctpgame.hpp/ctpgame.cpp`
- `ctpdata.hpp/ctpdata.cpp`
- `state.hpp/state.cpp`
- `strategy.hpp/strategy.cpp`
- `population.hpp/population.cpp`
- `moran.hpp/moran.cpp`
- `main.cpp`

# System Requirements

## Hardware Requirements

This code was built on a standard Macbook pro (Monterey 12.6) with 
- RAM: 8 GB
- CPU: 2.4 GHz Quad-Core Intel Core i5

## Software Requirements

OS Requirements

The package compilation was tested on OS X (laptop) and Unix (cluster) operating systems. 

To compile you will need
- GSL (tested for version 2.6 Unix platform and version 2.7 on mac os X) 
- GCC (version gcc-9.3) on Unix platform and clang on Mac os X.

# Compilation Guide
A CMAKE file is provide. Version 3.16.4 was used on the Unix platform. 
Compile by running first
```
cmake .
```
 to create the makefiles.  Once that works run 
```
make
```
which should produce an executable `ctpmoran`. 


# Demo
When executing `ctpmoran` without modifications in the `main.cpp`, the program will produce for an ICG with L=4 the results for β=0.3 and ε=0.18. It will run 100 times 10^7 iterations. One run for 10^8 iterations takes around 29000 seconds or 8 hours.

In `main.cpp` there are 13 parameters that can be set.These are
- `length` The length of the ICG (default value is 4).
- `first` The payoff player 1 gets in the first step (default value is 0.4).
- `second` The payoff player 2 gets in the first step (default value is 0.1).
- `factor` The growth of the resource at each step (default value is 2.0).
- `levels` Reasoning levels of the individuals (default value is 4).
- `maxlevel` Maximum reasoning level in the population (default value is 4).
- `epsilon` Reasonig error for each indivisual (default value is 0.18).
- `betas` Selection strength in the stochastic evolutionary dynamic (default value is 0.3).
- `repeats` Number of repeats to calculate the fitness of each individual in the monte carlo sampling process (default value is 50000).
- `psize` Population size (default value is 500).
- `mut` Mutation probability (default value is 0.0).
- `runs` Number of repititions of the simulation (default value is 100.0).
- `iterations` Number of steps of the moran process (default value is 10^7).
- `cost` Cost associated with each additional reasoning level (default value is 0.0).


In `main.cpp` there is 1 function that provides the results  to reproduce the data in the Extended Data Figure 8:
- `runMoranSimulations` : this function runs the moran simulation. 

# Output 
The raw output files as well as the R files needed to reproduce the figures can be found in the folder `results`.

The R-scripts can simply be excuted in the same folder as the raw data files.

# Citation

When using or extending this code, please cite :
T. Lenaerts, J.M. Pacheco and F.C. Santos (2022) Evolution of a Theory of Mind ...
See als Zenodo; https://zenodo.org/records/10229394 

