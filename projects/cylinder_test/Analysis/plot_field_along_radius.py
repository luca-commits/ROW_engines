import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import copy

mu = 1 * 4 * np.pi * 1e-7

# Load the data from the CSV file
# dataframe_004_spherical = pd.read_csv("transient/conductivity_10000.csv")
dataframe_004_spherical = pd.read_csv("transient/test1.csv")

# Filter the dataframe to include only rows where the arc_length is less than 10
dataframe_004_spherical = dataframe_004_spherical[dataframe_004_spherical['arc_length'] < 5]

current = 1
total_current = 1
# Initialize the analytic array with zeros, using the same index as the dataframe
analytic = copy.deepcopy(dataframe_004_spherical['arc_length']) * (mu * current / 2)

# Calculate the inverse of the arc_length for values where arc_length > circle_radius
circle_radius = 1

area = np.pi * circle_radius**2
current_density = total_current * area

arc_length = dataframe_004_spherical['arc_length']

# Create a boolean mask for the condition arc_length > circle_radius
mask = arc_length > circle_radius

# Calculate one_over_r for the arc_length values that meet the condition
analytic[mask] = (current * mu /2 ) / arc_length[mask] 

plt.Figure()
B_magnitude_004_spherical = (dataframe_004_spherical['B:0'].apply(np.square) + dataframe_004_spherical['B:1'].apply(np.square)).apply(np.sqrt)
B_magnitude_004_spherical = B_magnitude_004_spherical


plt.plot(arc_length, analytic, label="analytic solution")
plt.plot(dataframe_004_spherical['arc_length'], B_magnitude_004_spherical, label="numeric solution")
plt.legend()
plt.xlabel("radius along cylinder")
plt.ylabel("B field")
plt.savefig("B_field_for_various_meshes.png")
# plt.plot(dataframe['arc_length'], analytic_sq)
plt.show()
 