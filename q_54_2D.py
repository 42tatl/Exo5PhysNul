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

print("Parameters:")
print(f"cb_gauche: {cb_gauche}, cb_droite: {cb_droite}")
print(f"v_uniform: {v_uniform}, impose_nsteps: {impose_nsteps}")
print(f"A: {A}, f_hat: {f_hat}, x1: {x1}, x2: {x2}")
print(f"tfin: {tfin}, equation_type: {equation_type}")
print(f"nx: {nx}, n_init: {n_init}, initialization: {initialization}")
print(f"initial_state: {initial_state}, CFL: {CFL}")
print(f"nsteps: {nsteps}, output: {output}")
print(f"n_stride: {n_stride}, ecrire_f: {ecrire_f}")
print(f"hL: {hL}, hR: {hR}, h00: {h00}")
print(f"xa: {xa}, xb: {xb}, L: {L}, om: {om}")

cases = ["B"]
os.makedirs("outputs", exist_ok=True)

# Run simulations
for case in cases:
    output_name = f"wave_54_{case}"  
    fct.run_simulation(executable, input_filename, output_name, **params)

# Visualization
for case in cases:
    x, t, f, v, e = fct.read_wave_data(f"outputs/wave_54_{case}")
    # Check array orientation plots
    plt.figure(figsize=(8, 5))
    plt.imshow(f, aspect='auto',
                  extent=[t[0], t[-1], x[0], x[-1]],
                  origin='lower', cmap='turbo')
    plt.colorbar(label=r"$f(x,t)$ [m]")
    plt.xlabel(r"Time $t$ [s]")
    plt.ylabel(r"Position $x$ [m]")
    plt.title(f"Wave Propagation ({case})")
    plt.tight_layout()
    fct.save_figure(f"wave_{case}_54_x_vs_t.png")
    plt.show()


