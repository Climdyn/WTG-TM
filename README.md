# WTG-TM: The Weak Temperature Gradient forced Tropical Model

## About ##

(c) 2026 RMIB - Climdyn - Dynamical Meteorology and Climatology

See LICENSE for license information.

This software is provided as supplementary material with:

* Vannitsem, S., and Demaeyer, J.: Emergence of chaos in the tropical atmosphere: Study of the Weak Temperature 
Gradient system. Submitted to Nonlinear Processes in Geophysics.

**Please cite this article if you use (a part of) this software for a
publication.**

The authors would appreciate it if you could also send a reprint of
your paper to <jonathan.demaeyer@meteo.be> and
<svn@meteo.be>. 

Consult the WTG-TM [code repository](http://www.github.com/Climdyn/WTG-TM)
for updates, and [our website](http://climdyn.meteo.be) for additional
resources.

## Installation

The code requires at least Python 3.10 and the [LayerCake](https://github.com/Climdyn/LayerCake) framework.
In addition, the [qgs model](https://github.com/Climdyn/qgs) integrators and Lyapunov estimators are also used.
To install LayerCake and qgs, the simplest way is to use [pip](https://pypi.org/).
In a terminal, run:

```bash
pip install layercake-model
pip install qgs
```

and you are set.
In addition, the generated Fortran codes need to be compiled, for example with [GNU Fortran](https://gcc.gnu.org/fortran/).

## Codes description

### `generate_fortran_code.py`

Upon running

```bash
python generate_fortran_code.py
```
this code will generate a file `WTG-TM_Heun.f90` with by default the equations of the model truncated at wavenumber 2.
These equations are then integrated with the [Heun algorithm](https://en.wikipedia.org/wiki/Heun's_method) to give 
the time evolution of the model.
The Fortran code can then be compiled and run with:
```bash
gfortran WTG-TM_Heun.f90 && ./a.out
```
This result in a file `evol_field.dat` being produced, with the first column being the time, and the others being the 
spectral components $\psi_i$ of the fields. 
The default model resolution setup can be modified by editing the Python script. 
The parameters of the model can be changed by editing the generated Fortran code.
Finally, the Fortran compiler used here is `gfortran`, but other compilers should work as well.

### `generate_ode_latex_equation.py`

This code will generate and print the equations of the model in LaTeX.
It was used to get the equations appearing in the manuscript above.
By default, it will print the equations of the model truncated at wavenumber 2, but this can be changed by editing 
the script.

### `compute_trajectory.py`

This code computes a trajectory of the model and saves it like the Fortran code does, but in addition it rescales the 
time in days.
By default, it also generates a movie of the spatial fields of the model, saved in the file `movie.mp4`, but this can 
be deactivated.
Again, by default, it will compute the model truncated at wavenumber 2, but this can be changed by editing
the script. The forcing $\chi$  can also be changed using this script, while the more general model parameters can be
changed by editing the `model_definition.py` module.

### `compute_lyapunov_spectrum.py`

This script computes and plots the Lyapunov spectrum, using the qgs integrators.
By default, it will compute this for the model truncated at wavenumber 2, but this can be changed by editing
the script. The the forcing $\chi$ of the model can also be changed this way.
More general model parameters can be changed by editing the `model_definition.py` module.

