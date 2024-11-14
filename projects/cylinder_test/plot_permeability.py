import numpy as np
import matplotlib.pyplot as plt

# Read computed data
data = np.loadtxt('magnetic_data.txt', skiprows=1)
H = data[:, 0]
B = data[:, 1] * 1e6  # Convert from 1e-6 Tesla to Tesla
mu_r = data[:, 2] * 1e6  # Convert to actual relative permeability values

# Original data points
H_original = np.array([100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,2000])
mu_r_original = np.array([0.0004,0.0006,0.00130,0.00225,0.00225,0.0017,0.00159,0.00144,0.00131,0.00121,0.00104,0.00092,0.00083,0.00068]) * 1e6

# Calculate B for original points
mu_0 = 4 * np.pi * 1e-7  # Permeability of free space
B_original = mu_r_original * mu_0 * H_original

# Create figure with two y-axes
plt.figure(figsize=(12, 8))
fig, ax1 = plt.subplots(figsize=(12, 8))
ax2 = ax1.twinx()

# Plot styling
ax1.grid(True, color='gray', linestyle='-', alpha=0.2)
plt.title('B/H and μr/H Curves', pad=20, fontsize=14)

# Plot B/H curve and points (green, left axis)
line1 = ax1.plot(H, B, 'g-', linewidth=2, label='B/H Curve')
ax1.plot(H_original, B_original, 'go', markersize=8, label='B/H Original Points')

ax1.set_xlabel('Magnetising Force H (A/m)', fontsize=12)
ax1.set_ylabel('Flux Density B (Tesla)', fontsize=12)

# Plot μr/H curve and points (red, right axis)
line2 = ax2.plot(H, mu_r, 'r-', linewidth=2, label='μr/H Curve')
ax2.plot(H_original, mu_r_original, 'ro', markersize=8, label='μr/H Original Points')
ax2.set_ylabel('Relative Permeability μr', fontsize=12)

# Set axis limits
ax1.set_xlim(0, 2000)
ax1.set_ylim(0, 2.0)
ax2.set_ylim(0, 3000)

# Make ticks and labels more visible
ax1.tick_params(axis='both', which='major', labelsize=10)
ax2.tick_params(axis='both', which='major', labelsize=10)

# Combine legends from both axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=12)

# Print some statistics to verify the data
print(f"B range: {B.min():.2e} to {B.max():.2e} Tesla")
print(f"μr range: {mu_r.min():.2f} to {mu_r.max():.2f}")
print(f"H range: {H.min():.2f} to {H.max():.2f} A/m")
print("\nOriginal data ranges:")
print(f"B range: {B_original.min():.2e} to {B_original.max():.2e} Tesla")
print(f"μr range: {mu_r_original.min():.2f} to {mu_r_original.max():.2f}")
print(f"H range: {H_original.min():.2f} to {H_original.max():.2f} A/m")

# Adjust layout
plt.tight_layout()

# Save with high resolution
plt.savefig('magnetic_curves.png', dpi=300, bbox_inches='tight')
plt.show()