
##############################################################################################
#
# This script compute the Lyapunov spectrum of the Weak Temperature Gradient (WTG) tropical
# model using the LayerCake and qgs libraries
#
##############################################################################################

from model_definition import define_model

import numpy as np
import matplotlib.pyplot as plt

from qgs.integrators.integrator import RungeKuttaIntegrator
from qgs.toolbox.lyapunov import LyapunovsEstimator

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
    dt = 0.01
    integrator.integrate(0., 200000., dt, ic=ic, write_steps=0)
    time, ic = integrator.get_trajectories()
    print('Done.')

    # making the Lyapunov spectrum estimation
    lvint = LyapunovsEstimator()

    print('Computing the Lyapunov spectrum...')
    lvint.set_func(f, Df)
    lvint.compute_lyapunovs(0., 10000., 20000., 0.01, 0.01, ic, write_steps=10)
    btl, btraj, bexp, bvec = lvint.get_lyapunovs()
    print('Done.')

    # plotting the results

    # trajectories
    plt.figure(figsize=(10, 8))

    plt.plot(btl*27998/(24*3600), btraj[:5].T)
    plt.xlabel('time [day]')
    plt.ylabel('trajectories')

    plt.figure(figsize=(15, 4))

    mean_exp = np.mean(bexp, axis=-1)*(24*3600)/27998

    x_pos = np.arange(1., model_definition.ndim+1, 1)

    plt.plot(x_pos, mean_exp)

    plt.xticks(x_pos, map(str, range(1, model_definition.ndim+1, 1)))
    plt.axhline(ls='--', color='k', lw=0.75)

    plt.xlim(x_pos[0]-1., x_pos[-1]+1.)
    plt.ylim(np.min(mean_exp)-0.1, np.max(mean_exp)+0.1)

    plt.ylabel("Lyapunov exponent [day$^{-1}$]")
    plt.xlabel("Index of the Lyapunov exponent")

    plt.show()

