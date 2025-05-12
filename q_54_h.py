import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct
import matplotlib as mpl

mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 14,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

# QUESTION A
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
os.makedirs("outputs", exist_ok=True)

# Run simulations
for case in cases:
    output_name = f"wave_54_{case}"  
    fct.run_simulation(executable, input_filename, output_name, **params)

# Visualization
for case in cases:
    x, t, f, v, e = fct.read_wave_data(f"outputs/wave_54_{case}")
    h = np.full_like(x, hL)
    for i in range(len(x)):
        if 0 <= x[i] <= xa:
            h[i] = -hL
        elif xa < x[i] <= xb:
            h[i] = -(0.5*(hL + hR) + 0.5*(hL - hR)*np.cos(np.pi*(x[i] - xa)/(xb - xa)))
        elif xb < x[i] <= L:
            h[i] = -hR
plt.figure(figsize=(8, 5))
plt.plot(x, h, label=r"$h(x)$", color='blue')

# Add vertical dotted lines and right-shifted LaTeX labels
for xpos, label in zip([xa, xb, x1, x2], [r"$x_a$", r"$x_b$", r"$x_1$", r"$x_2$"]):
    plt.axvline(x=xpos, color='black', linestyle='dotted')
    plt.text(xpos + 0.01 * L, max(h) * 1.02, label, rotation=90,
             verticalalignment='bottom', horizontalalignment='left')

plt.xlabel(r"$x$ [m]")
plt.ylabel(r"Ocean depth $h$ [m]")
plt.tight_layout()
fct.save_figure(f"wave_{case}_54_h_vs_x.png")
plt.show()


