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

input_filename = "config_5_4.in"
params = fct.read_in_file(input_filename)
(cb_gauche, cb_droite, v_uniform, 
 A, f_hat, x1, x2, 
 tfin, equation_type,
 nx, n_init, initialization, initial_state,
 CFL, nsteps, impose_nsteps,
 output, n_stride, ecrire_f,
 hL, hR, h00, _, _, L, om) = fct.get_wave_params(params)

# Output cases to compare
output_cases = [
    ("wave_54_xa500_xb950", 500_000, 950_000),
    ("wave_54_xa600_xb950", 600_000, 950_000),
    ("wave_54_xa700_xb950", 700_000, 950_000),
    ("wave_54_xa800_xb950", 800_000, 950_000),
    ("wave_54_xa850_xb950", 850_000, 950_000),
    ("wave_54_xa900_xb950", 900_000, 950_000),
]

# Create figures once
fig1, crests_ax = plt.subplots(figsize=(8, 5))
fig2, amps_ax = plt.subplots(figsize=(8, 5))
fig3, vel_ax = plt.subplots(figsize=(8, 5))

x_th = np.linspace(0, L, 1000)
g = 9.81

for case, xa, xb in output_cases:
    x, t, f, v, e = fct.read_wave_data(f"outputs/{case}")
    if f.shape[0] == len(x):
        f = f.T  # shape should be (nt, nx)

    t_crete = []
    amplitudes = []
    x_selected = []

    for i in range(1, len(x) - 1):
        signal = f[:, i]
        imax = np.argmax(signal)
        if imax == 0 or imax == len(t) - 1:
            continue

        t_fit = t[imax-1:imax+2]
        f_fit = signal[imax-1:imax+2]
        a, b, c = np.polyfit(t_fit, f_fit, 2)
        t_peak = -b / (2 * a)
        f_peak = a * t_peak**2 + b * t_peak + c

        t_crete.append(t_peak)
        amplitudes.append(f_peak)
        x_selected.append(x[i])

    x_arr = np.array(x_selected)
    t_arr = np.array(t_crete)

    # Crest plot
    crests_ax.plot(x_arr, t_arr, label=case)

    amps_ax.plot(x_arr, amplitudes, label=case)

    # Velocity
    k = 30
    v_num = []
    x_mid = []
    for i in range(k, len(x_arr) - k):
        dx = x_arr[i + k] - x_arr[i - k]
        dt = t_arr[i + k] - t_arr[i - k]
        if dt == 0:
            v_num.append(np.nan)
        else:
            v_num.append(dx / dt)
            x_mid.append(x_arr[i])

    vel_ax.plot(x_mid, v_num, label=case)

# Final plot formatting
for ax, xlabel, ylabel, fname in [
    (crests_ax, r"$x_i$ [m]", r"$t_{\mathrm{cr\hat{e}te},i}$ [s]", "comparison_crests.png"),
    (amps_ax, r"$x_i$ [m]", r"Amplitude $f(x_i, t_{\mathrm{cr\hat{e}te},i})$ [m]", "comparison_amps.png"),
    (vel_ax, r"$x_i$ [m]", r"Vitesse $v(x)$ [m/s]", "comparison_vel.png"),
]:
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    fct.save_figure(fname)
    plt.show()
