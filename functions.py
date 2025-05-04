import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from matplotlib.animation import FuncAnimation

def read_in_file(filename):
    '''Reads in a file and returns the data as a dictionary with proper types'''
    variables = {}
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("#") and '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.split('!')[0].strip()  # Remove comments
                
                # Try to determine the type automatically
                if value.lower() in ('true', 'false'):
                    variables[key] = value.lower() == 'true'
                elif value.replace('.', '').replace('e-', '').replace('e+', '').isdigit():
                    if '.' in value or 'e' in value.lower():
                        variables[key] = float(value)
                    else:
                        variables[key] = int(value)
                else:
                    variables[key] = value
    return variables

def get_wave_params(params):
    '''Extracts all wave simulation parameters with proper typing'''
    # Boundary conditions (strings)
    cb_gauche = params.get("cb_gauche", "libre")
    cb_droite = params.get("cb_droite", "fixe")
    
    # Boolean parameters
    v_uniform = bool(params.get("v_uniform", True))
    impose_nsteps = bool(params.get("impose_nsteps", True))
    
    # Numerical parameters
    A = float(params.get("A", 0.1))
    f_hat = float(params.get("f_hat", 1.0))
    x1 = float(params.get("x1", 3.0))
    x2 = float(params.get("x2", 8.0))
    tfin = float(params.get("tfin", 1.0))
    nx = int(params.get("nx", 64))
    n_init = int(params.get("n_init", 3))
    CFL = float(params.get("CFL", 1.0))
    nsteps = int(params.get("nsteps", 40))
    
    # String parameters
    equation_type = str(params.get("equation_type", "B"))
    initialization = str(params.get("initialization", "mode"))
    initial_state = str(params.get("initial_state", "static"))
    output = str(params.get("output", "test.out"))
    
    # Other parameters
    n_stride = int(params.get("n_stride", 1))
    ecrire_f = int(params.get("ecrire_f", 1))
    hL = float(params.get("hL", 8000.0))
    hR = float(params.get("hR", 20.0))
    h00 = float(params.get("h00", 4.0))
    xa = float(params.get("xa", 450000.0))  # Using float instead of scientific notation
    xb = float(params.get("xb", 950000.0))
    L = float(params.get("L", 15.0))
    om = float(params.get("om", 0.1))
    
    return (
        cb_gauche, cb_droite, v_uniform,
        A, f_hat, x1, x2,
        tfin, equation_type,
        nx, n_init, initialization, initial_state,
        CFL, nsteps, impose_nsteps,
        output, n_stride, ecrire_f,
        hL, hR, h00, xa, xb, L, om
    )

def run_simulation(executable, input_filename, output_name, **params):
    '''Runs the simulation with the given parameters'''
    os.makedirs("outputs", exist_ok=True)
    
    # Build command
    cmd = [executable, input_filename]
    for key, value in params.items():
        if key != "output":
            cmd.append(f"{key}={value}")
    cmd.append(f"output=outputs/{output_name}")
    
    print("\nRunning command:", " ".join(cmd))
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Command executed successfully")
        
        # Check for files WITHOUT .out extension
        expected_files = [
            f"outputs/{output_name}_en",
            f"outputs/{output_name}_f",
            f"outputs/{output_name}_v",
            f"outputs/{output_name}_x"
        ]
        
        if all(os.path.exists(f) for f in expected_files):
            print("All output files created successfully")
            return True
        else:
            missing = [f for f in expected_files if not os.path.exists(f)]
            print(f"Missing files: {missing}")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error:\n{e.stderr}")
        return False

def read_wave_data(output_base):
    """
    Reads wave simulation output files (without .out extension)
    Returns: x, t, f, v, energy
    """
    try:
        # Read files WITHOUT .out extension
        x = np.loadtxt(f"{output_base}_x")
        v = np.loadtxt(f"{output_base}_v")
        f_data = np.loadtxt(f"{output_base}_f")
        energy = np.loadtxt(f"{output_base}_en")
        
        # Process wave height data
        t = []
        f = []
        for row in f_data:
            if len(row) == len(x) + 1:  # Time + spatial points
                t.append(row[0])
                f.append(row[1:])
        
        return np.array(x), np.array(t), np.array(f), np.array(v), np.array(energy)
        
    except Exception as e:
        print(f"Error reading data: {e}")
        return None, None, None, None, None

def save_figure(filename, fig=None, subfolder="figures", dpi=300, tight=True):
    """Saves a Matplotlib figure to a specified subfolder"""
    if fig is None:
        fig = plt.gcf()
    os.makedirs(subfolder, exist_ok=True)
    filepath = os.path.join(subfolder, filename)
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight' if tight else None)
    print(f"Figure saved to {filepath}")

def animate_wave_propagation(x, t, f, save_as=None):
    """Creates an animation of wave propagation"""
    fig, ax = plt.subplots(figsize=(10, 6))
    line, = ax.plot(x, f[0, :], 'b-')
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(np.min(f), np.max(f))
    ax.set_xlabel('Position x (m)')
    ax.set_ylabel('Wave height f(x,t) (m)')
    ax.set_title('Wave Propagation')
    ax.grid(True)
    
    def update(frame):
        line.set_ydata(f[frame, :])
        return line,
    
    ani = FuncAnimation(fig, update, frames=len(t), interval=50, blit=True)
    
    if save_as:
        try:
            # Try MP4 first, fall back to GIF if ffmpeg not available
            ani.save(save_as, writer='ffmpeg', fps=30)
        except:
            gif_path = os.path.splitext(save_as)[0] + '.gif'
            ani.save(gif_path, writer='pillow', fps=15)
            print(f"Saved animation as {gif_path} (MP4 unavailable)")
    
    plt.show()
    return ani

def calculate_analytical_mode(x, n, L, t, u):
    """Calculates analytical solution for normal modes"""
    k = n * np.pi / L
    omega = u * k
    return np.sin(k * x) * np.cos(omega * t)