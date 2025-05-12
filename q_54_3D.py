# plot_wave_surface.py

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plots
from matplotlib import cm
import os
import functions as fct

# Set LaTeX-style font globally
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 11,
})

# File setup
executable = './Exe1.exe'
repertoire = r"C:\Users\Avril\Desktop\Exo5PhysNul"
os.chdir(repertoire)

input_filename = "config_5_4.in"
params = fct.read_in_file(input_filename)
(cb_gauche, cb_droite, v_uniform, 
 A, f_hat, x1, x2, 
 tfin, equation_type,
 nx, n_init, initialization, initial_state,
 CFL, nsteps, impose_nsteps,
 output, n_stride, ecrire_f,
 hL, hR, h00, xa, xb, L, om) = fct.get_wave_params(params)

cases = ["B"]

# Create output directory
os.makedirs("outputs", exist_ok=True)

# Run simulations
for case in cases:
    output_name = f"wave_54_{case}"
    fct.run_simulation(executable, input_filename, output_name, **params)

# Generate 3D surface plot
for case in cases:
    x, t, f, v, e = fct.read_wave_data(f"outputs/wave_54_{case}")

    # Transpose f to match the (x, t) meshgrid layout
    f = f.T  # Now shape is (69, 10001)

    T, X = np.meshgrid(t, x)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Surface plot
    surf = ax.plot_surface(X, T, f, cmap='turbo', edgecolor='none', antialiased=True)

    # Labels
    ax.set_xlabel(r"$x$ [m]")
    ax.set_ylabel(r"$t$ [s]")
    ax.set_zlabel(r"$f(x,t)$ [m]")

    # Colorbar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.7, pad=0.1)
    cbar.set_label(r"$f(x,t)$ [m]")

    # Optional view angle
    ax.view_init(elev=30, azim=-120)

    # Subfigure label
    ax.text2D(0.5, -0.12, r"\textbf{(b)}", transform=ax.transAxes,
              ha='center', va='center', fontsize=14)

    fct.save_figure(f"wave_{case}_54_surface_plot.png")
    plt.show()


