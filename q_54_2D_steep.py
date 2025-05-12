import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct

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

x_pairs = [(500_000, 950_000), (600_000, 950_000), (700_000, 950_000), (800_000, 950_000), (850_000, 950_000), (900_000, 950_000)]

os.makedirs("outputs", exist_ok=True)

for xa_val, xb_val in x_pairs:
    output_name = f"wave_54_xa{xa_val//1000}_xb{xb_val//1000}"
    params['xa'] = xa_val
    params['xb'] = xb_val
    fct.run_simulation(executable, input_filename, output_name, **params)

# Visualization
for xa_val, xb_val in x_pairs:
    x, t, f, v, e = fct.read_wave_data(f"outputs/wave_54_xa{xa_val//1000}_xb{xb_val//1000}")
    # Check array orientation plots
    plt.figure(figsize=(8, 5))
    plt.imshow(f, aspect='auto',
                  extent=[t[0], t[-1], x[0], x[-1]],
                  origin='lower', cmap='turbo')
    plt.colorbar(label=r"$f(x,t)$ [m]")
    plt.xlabel(r"Time $t$ [s]")
    plt.ylabel(r"Position $x$ [m]")
    plt.tight_layout()
    fct.save_figure(f"wave_54_xa{xa_val//1000}_xb{xb_val//1000}.png")
    plt.show()


