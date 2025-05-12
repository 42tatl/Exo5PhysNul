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
 hL, hR, h00, xa, xb, L, om) = fct.get_wave_params(params)

# Use pre-generated output files

output_cases = [
    ("wave_54_xa500_xb950", 500_000, 950_000),
    ("wave_54_xa600_xb950", 600_000, 950_000),
    ("wave_54_xa700_xb950", 700_000, 950_000),
    ("wave_54_xa800_xb950", 800_000, 950_000),
    ("wave_54_xa850_xb950", 850_000, 950_000),
    ("wave_54_xa900_xb950", 900_000, 950_000),
]

for case, xa, xb in output_cases:
    x, t, f, v, e = fct.read_wave_data(f"outputs/{case}")
    # Transpose f to match the (x, t) meshgrid layout
    # f shape is (nt, nx) or (nx, nt), adjust accordingly

    if f.shape[0] == len(x):
        f = f.T  # ensure shape is (nt, nx)

    t_crete = []
    amplitudes = []
    x_selected = []

    for i in range(1, len(x) - 1):
        signal = f[:, i]  # f(t, x_i)
        
        # Find index of maximum
        imax = np.argmax(signal)
        
        # Skip boundaries
        if imax == 0 or imax == len(t) - 1:
            continue
        
        # Quadratic fit: use points around the maximum
        t_fit = t[imax-1:imax+2]
        f_fit = signal[imax-1:imax+2]
        
        # Fit quadratic: f(t) ≈ a*t^2 + b*t + c
        coeffs = np.polyfit(t_fit, f_fit, 2)
        a, b, c = coeffs
        
        # Vertex of parabola: t_crete = -b/(2a)
        t_peak = -b / (2 * a)
        f_peak = a * t_peak**2 + b * t_peak + c
        
        # Store values
        t_crete.append(t_peak)
        amplitudes.append(f_peak)
        x_selected.append(x[i])

    #CRESTS
    plt.figure(figsize=(8, 5))
    plt.plot(x_selected, t_crete)
    plt.xlabel(r"$x_i$ [m]")
    plt.ylabel(r"$t_{\mathrm{cr\hat{e}te},i}$ [s]")
    plt.grid(True)
    plt.tight_layout()
    fct.save_figure(f"{case}_crests.png")
    plt.show()

    #AMPLITUDES
    #WKB
    x_th = np.linspace(0, L, 1000)
    h_th = np.full_like(x_th, hL)

    for i in range(len(x_th)):
        if 0 <= x_th[i] <= xa:
            h_th[i] = hL
        elif xa < x_th[i] <= xb:
            h_th[i] = 0.5*(hL + hR) + 0.5*(hL - hR)*np.cos(np.pi*(x_th[i] - xa)/(xb - xa))
        elif xb < x_th[i] <= L:
            h_th[i] = hR

    A_th = np.sqrt(np.sqrt(hL / h_th))
    plt.figure(figsize=(8, 5))
    plt.plot(x_selected, amplitudes)
    plt.plot(x_th, A_th, label=r"$A_{\mathrm{WKB}}(x)$", color='red', linestyle='dashed')
    plt.xlabel(r"$x_i$ [m]")
    plt.ylabel(r"Amplitude $f(x_i, t_{\mathrm{cr\hat{e}te},i})$ [m]")
    plt.grid(True)
    plt.tight_layout()
    fct.save_figure(f"{case}_amps.png")
    plt.show()

    #Velocity calculation
    k = 30
    x_arr = np.array(x_selected)
    t_arr = np.array(t_crete)

    v = []
    x_mid = []

    for i in range(k, len(x_arr) - k):
        dx = x_arr[i + k] - x_arr[i - k]
        dt = t_arr[i + k] - t_arr[i - k]
        if dt == 0:
            v.append(np.nan)
        else:
            v.append(dx / dt)
            x_mid.append(x_arr[i])  # center point

    #WKB
    g=9.81  # m/s^2
    v_th = np.sqrt(g * h_th)  # Assuming g = 9.81 m/s^2

    plt.plot(x_mid, v)
    plt.plot(x_th, v_th, label=r"$v_{\mathrm{WKB}}(x)$", color='red', linestyle='dashed')
    plt.xlabel(r"$x_i$ [m]")
    plt.ylabel(r"Vitesse $v(x)$ [m/s]")
    plt.grid(True)
    plt.tight_layout()
    fct.save_figure(f"{case}_vel.png")
    plt.show()
