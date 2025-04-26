import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from scipy.integrate import simpson

def analyze_power_loss_error(file1_path, file2_path, time_step=1.0):
    """
    Analyze the error between two power loss simulation files.
    Computes the sum of squared differences and integrates the error over time.

    Args:
        file1_path (str): Path to the first CSV file
        file2_path (str): Path to the second CSV file
        time_step (float): Time step between consecutive rows (default=1.0)

    Returns:
        tuple: (float: Sum of squared differences, float: Integrated error over time)
    """
    try:

        data1 = pd.read_csv(file1_path, header=None, skiprows=1, delimiter=',')
        data2 = pd.read_csv(file2_path, header=None, skiprows=1, delimiter=',')


        data1 = data1.fillna(0).astype(float)
        data2 = data2.fillna(0).astype(float)

        array1 = data1.values
        array2 = data2.values


        num_steps1 = array1.shape[0]
        num_steps2 = array2.shape[0]

        min_steps = min(num_steps1, num_steps2)
        array1 = array1[:min_steps]
        array2 = array2[:min_steps]


        if array1.shape != array2.shape:
            raise ValueError(f"Files have different shapes after truncation: {array1.shape} vs {array2.shape}")
    

        difference = np.abs(array2 - array1)
        norm = np.abs(np.sum(array1, axis=1))
        difference = np.sum(difference, axis=1)
        relative_error = difference / norm
        print("Relative Error", relative_error)
        error_over_time = np.mean(relative_error, axis=0)
        print(f"Average error over time: {error_over_time}")


    except Exception as e:
        print(f"Error: {e}")
        return None, None

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py file1.csv file2.csv")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    analyze_power_loss_error(file1_path, file2_path)
