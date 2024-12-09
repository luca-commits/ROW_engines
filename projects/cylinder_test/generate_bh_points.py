import numpy as np
import matplotlib.pyplot as plt

def sample_magnetic_field(mu_0=1.256637e-6, num_points=1000):
    k1 = 4.847
    k2 = 1.908
    k3 = 227.3
    # Create arrays starting from a small positive value
    x = np.logspace(-7, 0, num_points)
    # f_x = mu_0 * (np.arctan(x * 100) + x)
    f_x = ( 1/(k1 np.exp(k2 * x * x) + k3)) + mu_0
    
    # Add (0,0) point by prepending it to the arrays
    x = np.concatenate(([0], x))
    f_x = np.concatenate(([0], f_x))
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.scatter(f_x, x, s=1, label='B-H Curve')
    plt.ylabel('H (Magnetic Field Strength)')
    plt.xlabel('B (Magnetic Flux Density)')
    plt.title('B-H Curve with Enhanced Initial Region Sampling')
    plt.grid(True)
    plt.legend()
    
    # Set axis limits
    # plt.xlim(-0.1, 1)
    # plt.ylim(-0.1, 1)
    plt.show()
    
    # Save arrays to file in the requested format
    with open('tabular_data.txt', 'w') as f:
        f.write(f"[{','.join(f'{val:.20f}' for val in f_x)}]\n")
        f.write(f"[{','.join(f'{val:.20f}' for val in x)}]")
    
    return f_x, x

if __name__ == "__main__":
    # Generate and save samples
    B_values, H_values = sample_magnetic_field()
    print("Data has been saved to 'magnetic_field_data.txt'")