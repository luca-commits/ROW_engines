import numpy as np

def get_reluctivity(B, max_=5000, c_=100, mu0_=4 * np.pi * 1e-7):
    if B < 1e-12:
        B = 1e-12
    relative_permeability = max_ / (1 + np.power(B, 4) * max_ / c_) + 1
    return (1 / relative_permeability) / mu0_

# Generate B values (0 to 2 Tesla with small steps)
B_values = np.linspace(0, 2, 100)
H_values = np.array([B * get_reluctivity(B) for B in B_values])

# Save to CSV
np.savetxt('hb_curve.csv', np.column_stack((H_values, B_values)), delimiter=',', header='H,B', comments='')
