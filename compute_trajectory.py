
##############################################################################################
#
# This script compute a trajectory of the Weak Temperature Gradient (WTG) tropical
# model using the LayerCake and qgs libraries.
# Optionally it can also generate a movie.
#
##############################################################################################

from model_definition import define_model

import numpy as np
import matplotlib.pyplot as plt

from qgs.integrators.integrator import RungeKuttaIntegrator

video = True

# Guarding the main script to deal with multiprocessing import issue
# in case the start method is 'spawn' or 'forkserver'.
# See https://docs.python.org/3/library/multiprocessing.html#the-spawn-and-forkserver-start-methods
# 'Safe importing of main module' section for more details.
if __name__ == "__main__":

    print('Creating the model...')
    # constructing the model
    model_definition = define_model(2, 2)

    # computing the tensor (might take a long time depending on the resolution
    model_definition.compute_tensor(numerical=True, compute_inner_products=True)

    # generating the tendencies and Jacobian Numba callable
    f, Df = model_definition.compute_tendencies()
    print('Done.')

    print('Integrating to get an initial condition on the attractor...')
    # defining a RK4 integrator (from qgs)
    integrator = RungeKuttaIntegrator()
    integrator.set_func(f)

    # integrating to get a first initial condition on the attractor
    ic = np.random.rand(model_definition.ndim)*0.001
    dt = 0.01  # timestep
    integrator.integrate(0., 200000., dt, ic=ic, write_steps=0)
    _, ic = integrator.get_trajectories()
    print('Done.')

    print('Integrating to get a trajectory on the attractor...')
    # integrating to get a first initial condition on the attractor
    integrator.integrate(0., 20000., dt, ic=ic, write_steps=10)
    time, trajectory = integrator.get_trajectories()
    data = np.concatenate((27998 * time[np.newaxis, ...] / (24 * 3600), trajectory))
    np.savetxt('WTG_TM_trajectory.dat', data.T)
    print('Done.')

    if video:

        # Domain specification
        # Geometry / modes
        n = 0.20
        xmax = 2 * np.pi / n

        # Grid
        nx = 31
        ny = 21
        x = np.linspace(0, xmax, nx)
        y = np.linspace(-np.pi / 2, np.pi / 2, ny)
        X, Y = np.meshgrid(x, y)

        # Basis function specification
        basis = model_definition.layers[0].equations[0].terms[2].field.basis
        blist = basis.num_functions()
        dblist = basis.directional_derivative()





