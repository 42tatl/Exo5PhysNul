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

'''
# Run simulations
for case in cases:
    current_params = params.copy()
    current_params["initial_state"] = case
    output_name = f"wave_{case}"
    fct.run_simulation(executable, "config_5_3a.in", output_name, **current_params)


for case in cases:
    x, t, f, _, _ = fct.read_wave_data(f"outputs/wave_" + case)
    if x is not None:
        f = np.array(f)
        # Vérification de l’ordre (on veut f[x, t])
        if f.shape[0] == len(t):  # alors f = f[t, x], il faut transposer
            f_plot = f.T
        else:
            f_plot = f

        plt.figure(figsize=(8, 5))
        plt.imshow(f_plot, aspect='auto',
                   extent=[t[0], t[-1], x[0], x[-1]],
                   origin='lower', cmap='turbo')
        plt.colorbar(label=r"$f(x,t)$ [m]")
        plt.xlabel(r"Time $t$ [s]")
        plt.ylabel(r"Position $x$ [m]")
        plt.tight_layout()
        fct.save_figure(f"wave_{case}_x_vs_t.png")
        plt.show()


for case in cases:
    x, t, f, _, _ = fct.read_wave_data(f"outputs/wave_" + case)
    if x is not None:
        f = np.array(f)
        
        # Vérification de l’ordre (on veut f[x, t])
        if f.shape[0] == len(t):  # alors f = f[t, x], il faut transposer
            f_plot = f.T
        else:
            f_plot = f
        plt.figure(figsize=(8, 5))
        plt.imshow(f_plot.T, aspect='auto',  # transpose pour avoir t en y
                   extent=[x[0], x[-1], t[0], t[-1]],
                   origin='lower', cmap='turbo')
        plt.colorbar(label=r"$f(x,t)$ [m]")
        plt.xlabel(r"Position $x$ [m]")
        plt.ylabel(r"Time $t$ [s]")
        plt.tight_layout()
        fct.save_figure(f"wave_{case}_t_vs_x.png")
        plt.show()
'''




#Question b for left initial state


# Valeurs de CFL à tester
beta_values = [0.5, 1.0, 1.0003]  # On ne teste que les valeurs stables
cases = ["left"]  

os.makedirs("outputs", exist_ok=True)

# Boucle sur les valeurs de β_CFL
for beta in beta_values:
    current_params = params.copy()
    current_params["CFL"] = beta
    output_name = f"wave_beta_{beta}"
    fct.run_simulation(executable, "config_5_3a.in", output_name, **current_params)

# Visualisation des résultats
plt.figure(figsize=(12, 8))

for beta in beta_values:
    x, t, f, _, _ = fct.read_wave_data(f"outputs/wave_beta_{beta}")
    if x is not None:
        f = np.array(f)
        if f.shape[0] == len(t):  # f est (t, x), on veut (x, t)
            f_plot = f.T
        else:
            f_plot = f

        # On affiche un snapshot à un temps fixe (par exemple t = tfin/2)
        t_idx = len(t) // 2
        plt.plot(x, f_plot[:, t_idx], label=f"β = {beta}")

plt.xlabel(r"Position $x$ [m]")
plt.ylabel(r"$f(x, t=t_{fin}/2)$ [m]")
plt.legend()
plt.grid(True)
plt.tight_layout()
fct.save_figure("stability_limit_beta.png")
plt.show()




