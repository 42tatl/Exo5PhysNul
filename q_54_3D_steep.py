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


output_cases = [
    "wave_54_xa500_xb950",
    "wave_54_xa600_xb950",
    "wave_54_xa700_xb950",
    "wave_54_xa800_xb950",
    "wave_54_xa850_xb950",
    "wave_54_xa900_xb950",
]

# Generate 3D surface plots
for case in output_cases:
    x, t, f, v, e = fct.read_wave_data(f"outputs/{case}")

    f = f.T  # Shape: (space, time)
    T, X = np.meshgrid(t, x)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, T, f, cmap='turbo', edgecolor='none', antialiased=True)

    ax.set_xlabel(r"$x$ [m]")
    ax.set_ylabel(r"$t$ [s]")
    ax.set_zlabel(r"$f(x,t)$ [m]")

    cbar = fig.colorbar(surf, ax=ax, shrink=0.7, pad=0.1)
    cbar.set_label(r"$f(x,t)$ [m]")

    ax.view_init(elev=30, azim=-120)

    ax.text2D(0.5, -0.12, r"\textbf{(b)}", transform=ax.transAxes,
              ha='center', va='center', fontsize=14)

    fct.save_figure(f"{case}_surface_plot.png")
    plt.show()
