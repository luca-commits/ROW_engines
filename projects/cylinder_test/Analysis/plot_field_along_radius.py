import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import copy

dataframe = pd.read_csv("B_field_along_r.csv")

circle_radius = 1


analytic = copy.deepcopy(dataframe['arc_length']) * 5e6
one_over_r = 2e13 / (analytic )
analytic[dataframe['arc_length'] > circle_radius] = one_over_r[dataframe['arc_length'] > circle_radius]

analytic_sq = copy.deepcopy(dataframe['arc_length']) * 5e6
one_over_r_sq = 2e19 / (analytic_sq )**2
analytic_sq[dataframe['arc_length'] > circle_radius] = one_over_r_sq[dataframe['arc_length'] > circle_radius]


plt.Figure()
B_magnitude = (dataframe['B:0'].apply(np.square) + dataframe['B:1'].apply(np.square)).apply(np.sqrt)
plt.plot(dataframe['arc_length'], B_magnitude)
plt.plot(dataframe['arc_length'], analytic)
plt.plot(dataframe['arc_length'], analytic_sq)
plt.show()
 