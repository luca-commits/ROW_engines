import numpy as np
import matplotlib.pyplot as plt

def read_data_points(filename='magnetic_field_data.txt'):
    with open(filename, 'r') as f:
        lines = f.readlines()
        B_points = np.array([float(x) for x in lines[0].strip()[1:-1].split(',')])
        H_points = np.array([float(x) for x in lines[1].strip()[1:-1].split(',')])
    return B_points, H_points

def plot_magnetic_data(data_file='magnetic_data.txt', points_file='magnetic_field_data.txt'):
    # Read interpolation points for reference
    B_points, H_points = read_data_points(points_file)
    
    # Read main interpolated data
    data = np.loadtxt(data_file, skiprows=1)
    B, H, nu, dnu = data.T
    
    # Calculate permeability (μ = 1/ν)
    mu = 1/nu
    
    # Find first non-nan index
    start_idx = np.where(~np.isnan(nu))[0][0]
    
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(24, 6))
    
    # B vs H
    ax1.plot(B, H, 'b-', linewidth=2, label='Interpolated')
    ax1.scatter(B_points, H_points, color='red', s=50, label='Reference points')
    ax1.set_xlabel('B [T]')
    ax1.set_ylabel('H [A/m]')
    ax1.set_title('B vs H')
    ax1.grid(True)
    ax1.legend()
    ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    
    # Reluctivity
    ax2.plot(B[start_idx:], nu[start_idx:], 'r-', linewidth=2)
    ax2.set_xlabel('B [T]')
    ax2.set_ylabel('ν [m/H]')
    ax2.set_title('Reluctivity vs B')
    ax2.grid(True)
    ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    # Derivative
    ax3.plot(B[start_idx:], dnu[start_idx:], 'g-', linewidth=2)
    ax3.set_xlabel('B [T]')
    ax3.set_ylabel('dν/dB [m/H/T]')
    ax3.set_title('Reluctivity Derivative vs B')
    ax3.grid(True)
    ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    
    # Permeability
    ax4.plot(B[start_idx:], mu[start_idx:], 'm-', linewidth=2)
    ax4.set_xlabel('B [T]')
    ax4.set_ylabel('μ [H/m]')
    ax4.set_title('Permeability vs B')
    ax4.grid(True)
    ax4.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    
    plt.tight_layout()
    plt.show()

plot_magnetic_data()