import numpy as np
import matplotlib.pyplot as plt
import os
import functions as fct

#Question a
executable = './Exe'  # Remove .exe for Mac
repertoire = r"/Users/lilimouelle/Desktop/PHYSNUM/Exo5PhysNul"  # Modify for correct directory
output_template = "wave_{direction}.out"
os.chdir(repertoire)

input_filename = "config_5_3a.in"
params = fct.read_in_file(input_filename)
(cb_gauche, cb_droite, v_uniform, 
 A, f_hat, x1, x2, 
 tfin, equation_type,
 nx, n_init, initialization, initial_state,
 CFL, nsteps, impose_nsteps,
 output, n_stride, ecrire_f,
 hL, hR, h00, xa, xb, L, om) = fct.get_wave_params(params)

#cases = ["left", "right", "static"]
cases = ["left"]


os.makedirs("outputs", exist_ok=True)

# Run simulations
for case in cases:
    current_params = params.copy()
    current_params["initial_state"] = case
    output_name = f"wave_{case}"
    fct.run_simulation(executable, "config_5_3a.in", output_name, **current_params)



# Affichage spatio-temporel : x en ordonnée, t en abscisse
for case in cases:
    x, t, f, _, _ = fct.read_wave_data(f"outputs/wave_{case}")
    if x is not None:
        plt.figure(figsize=(8, 5))
        # imshow attend f sous forme (espace, temps), donc f déjà au bon format
        plt.imshow(f, aspect='auto',
                   extent=[t[0], t[-1], x[0], x[-1]],
                   origin='lower', cmap='turbo')
        plt.colorbar(label=r"$f(x,t)$ [m]")
        plt.xlabel(r"Time $t$ [s]")
        plt.ylabel(r"Position $x$ [m]")
        plt.tight_layout()
        fct.save_figure(f"wave_{case}_spacetime_x_t.png")
        plt.show()


#Question b for left initial state



