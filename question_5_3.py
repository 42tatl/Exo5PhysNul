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

cases = ["left", "right", "static"]
os.makedirs("outputs", exist_ok=True)

# Run simulations
for case in cases:
    current_params = params.copy()
    current_params["initial_state"] = case
    output_name = f"wave_{case}"
    fct.run_simulation(executable, "config_5_3a.in", output_name, **current_params)

# Plot results
plt.figure(figsize=(10,6))
for case in cases:
    # Read files WITHOUT .out extension
    x, t, f, _, _ = fct.read_wave_data(f"outputs/wave_{case}")
    if x is not None:
        plt.plot(x, f[0,:], label=f"Initial {case}")
        plt.plot(x, f[-1,:], '--', label=f"Final {case}")

plt.xlabel('Position x (m)')
plt.ylabel('Wave height f(x,t) (m)')
plt.title('Wave Propagation Results')
plt.legend()
plt.grid(True)
plt.savefig("wave_results.png")
plt.show()



