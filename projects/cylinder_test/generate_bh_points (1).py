import numpy as np
import matplotlib.pyplot as plt

def sample_magnetic_field(mu_0=1.256637e-6, num_points=1000):
    x = np.linspace(0, 200, num_points)
    
    # Scale the arctan to match the input range approximately
    # f_x = 5e-13 * np.arctan(x/1e-13) / (np.pi/2) + mu_0 * x # Scaled arctan
    f_x =  60 * np.arctan(x / 60) / (np.pi/2) + mu_0 * x # Scaled arctan
    
    # Create the plot with similar ranges
    plt.figure(figsize=(10, 6))
    plt.scatter(f_x, x, s=1, label='B-H Curve')
    plt.ylabel('H (Magnetic Field Strength)')
    plt.xlabel('B (Magnetic Flux Density)')
    plt.title('B-H Curve with Enhanced Initial Region Sampling')
    plt.grid(True)
    plt.legend()
    
    # Set similar axis limits
    plt.xlim(0, 100)
    plt.ylim(0, 300)
    plt.show()
    # Save arrays to file in the requested format
    with open('magnetic_field_data.txt', 'w') as f:
        f.write(f"[{','.join(f'{val:.20f}' for val in f_x)}]\n")
        f.write(f"[{','.join(f'{val:.20f}' for val in x)}]")
    
    return f_x, x

if __name__ == "__main__":
    # Generate and save samples
    B_values, H_values = sample_magnetic_field()
    print("Data has been saved to 'magnetic_field_data.txt'")