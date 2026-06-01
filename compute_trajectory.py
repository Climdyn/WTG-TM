
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
from matplotlib.animation import FuncAnimation, FFMpegWriter

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
        # Generating a movie of the evolving spatial fields
        # from the trajectory in the spectral space
        print("Generating a movie of the evolving spatial fields...")
        # Use 12 snapshots every 15 steps
        step = 5
        nframes = 300

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
        basis_funcs_list = basis.num_functions()
        m = np.zeros((model_definition.ndim, *X.shape))
        for j, func in enumerate(basis_funcs_list):
            m[j] = func(X, Y)

        # Derivatives of basis functions
        dbasis = basis.directional_derivative()
        basis_dx_funcs_list = dbasis['x'].num_functions()
        dm_dx = np.zeros((model_definition.ndim, *X.shape))
        for j, func in enumerate(basis_dx_funcs_list):
            dm_dx[j] = func(X, Y)

        basis_dy_funcs_list = dbasis['y'].num_functions()
        dm_dy = np.zeros((model_definition.ndim, *X.shape))
        for j, func in enumerate(basis_dy_funcs_list):
            dm_dy[j] = func(X, Y)

        # Projecting the forcing Chi onto the spatial space
        chi_field = model_definition.layers[0].equations[0].terms[-2].terms[1].field
        chi_i = np.array(chi_field.parameters, dtype=float)
        chi = np.tensordot(chi_i, m, axes=1)
        dchi_dx = np.tensordot(chi_i, dm_dx, axes=1)
        dchi_dy = np.tensordot(chi_i, dm_dy, axes=1)

        # Projecting the wind and the streamfunction onto the spatial space
        np.array(chi_field.parameters, dtype=float)
        psis, Us, Vs = [], [], []
        tvals = []
        for j, state in enumerate(trajectory.T[:nframes*step:step]):
            psi = np.tensordot(state, m, axes=1)
            dpsi_dx = np.tensordot(state, dm_dx, axes=1)
            dpsi_dy = np.tensordot(state, dm_dy, axes=1)

            u = -dpsi_dy + dchi_dx
            v = dpsi_dx + dchi_dy

            psis.append(psi)
            Us.append(u)
            Vs.append(v)
            tvals.append(27998 * time[j*step] / (24 * 3600))

        # Generating the video
        vmin = np.min(psis)
        vmax = np.max(psis)
        levels = [-0.1, -0.05, -0.03, -0.01, 0.01, 0.03, 0.05, 0.1]

        fig, ax = plt.subplots(figsize=(8.8, 4.8))
        im = ax.imshow(
            psis[0],
            aspect="auto",
            origin="lower",
            extent=(0., xmax, y[0], y[-1]),
            vmin=vmin,
            vmax=vmax,
        )
        ax.axhline(0, color="black", linewidth=0.8)
        ax.contour(X, Y, chi, levels=levels, colors="white", linewidths=0.8)
        quiv = ax.quiver(X, Y, Us[0], Vs[0], color="white", scale=None, width=0.0025)

        ax.set_xticks([0, xmax / 2, xmax])
        ax.set_xticklabels(["0", r"$\pi/n$", r"$2\pi/n$"])
        ax.set_yticks([-np.pi / 2, 0, np.pi / 2])
        ax.set_yticklabels([r"$-\pi/2$", "0", r"$\pi/2$"])
        title = ax.set_title(rf"Velocity from $\psi$ and $\Xi$  (t = {tvals[0]:.3f} days)")

        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(r"$\psi$")


        def update(i):
            im.set_data(psis[i])
            quiv.set_UVC(Us[i], Vs[i])
            title.set_text(rf"Wind and streamfunction at t=({tvals[i]:.3f} days)")
            return [im, quiv, title]


        anim = FuncAnimation(fig, update, frames=len(psis), interval=500, blit=False)

        mp4_path = "movie.mp4"
        writer = FFMpegWriter(fps=4, bitrate=1800)
        anim.save(mp4_path, writer=writer)
        plt.close(fig)

        print("Done.")
        print("Saved animation:", mp4_path)
        print("Frames:", len(psis))
        print("Times used (in days):", [float(v) for v in tvals])







