# Re-import libraries after state reset

import numpy as np

import matplotlib.pyplot as plt



# Define parameters

B_saturation = 1.0   # Saturation flux density

mu_initial = 0.5     # Initial permeability

alpha = 0.8          # Scaling factor



# Define the function

def magnetic_flux_density(H_magnitude, B_saturation, mu_initial, alpha):

    return B_saturation * np.tanh(mu_initial * H_magnitude / (B_saturation * alpha))



# Generate H_magnitude values (only for x > 0)

H_magnitude = np.linspace(0, 2, 500)

B_curve = magnetic_flux_density(H_magnitude, B_saturation, mu_initial, alpha)



# Plot the curve

plt.figure(figsize=(8, 5))

plt.plot(H_magnitude, B_curve, label=f"B_saturation={B_saturation}, mu_initial={mu_initial}, alpha={alpha}", color="blue")

plt.axhline(B_saturation, color="red", linestyle="--", label="Saturation Limit")

plt.axvline(0, color="black", linestyle=":", linewidth=0.7)

plt.axhline(0, color="black", linestyle=":", linewidth=0.7)

plt.title("Magnetic Flux Density Curve (x > 0)", fontsize=14)

plt.xlabel("H_magnitude (Magnetic Field Intensity)", fontsize=12)

plt.ylabel("B (Magnetic Flux Density)", fontsize=12)

plt.legend(fontsize=10)

plt.grid(alpha=0.5)

plt.show()
