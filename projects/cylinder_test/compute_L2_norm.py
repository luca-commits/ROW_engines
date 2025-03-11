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
        # Read the data from the CSV files, skipping the first line (comment)
        data1 = pd.read_csv(file1_path, header=None, skiprows=1, delimiter=',')
        data2 = pd.read_csv(file2_path, header=None, skiprows=1, delimiter=',')

        print(data1)

        # Fill NaN values with 0 and convert to numerical types
        # data1 = data1.fillna(0).astype(float)
        # data2 = data2.fillna(0).astype(float)

        # Convert to numpy arrays for easier computation
        array1 = data1.values
        array2 = data2.values

        # Determine the number of time steps in each file
        num_steps1 = array1.shape[0]
        num_steps2 = array2.shape[0]

        # Truncate the arrays to the minimum number of time steps
        min_steps = min(num_steps1, num_steps2)
        array1 = array1[:min_steps]
        array2 = array2[:min_steps]

        # Check if the arrays have the same shape after truncation
        if array1.shape != array2.shape:
            raise ValueError(f"Files have different shapes after truncation: {array1.shape} vs {array2.shape}")

        print("Minimum value of abs(array2):", np.min(np.abs(array2)))
        print("Number of zeros in abs(array2):", np.sum(np.abs(array2) == 0))
        
        tolerance = 1e-10  # Adjust as needed

        squared_diff = np.divide(np.abs(array1 - array2), np.abs(array2), out=np.zeros_like(array2), where=array2 > tolerance)



        np.set_printoptions(threshold=np.inf, linewidth=np.inf)  # Print full arrays

        print("difference" , (array1 - array2)[:10])
        print("rhs" , array2[:10])

        # Compute the sum of all squared differences
        sum_of_squares = np.mean(squared_diff)

        # Calculate error at each time point (assuming each row is a time point)
        error_over_time = np.mean(squared_diff, axis=1)  # Sum across all spatial points for each time

        # Create time array based on the number of rows
        time_points = np.arange(0, len(error_over_time) * time_step, time_step)

        # Integrate error over time using Simpson's rule
        if len(error_over_time) > 1:
            integrated_error = simpson(error_over_time, time_points)
        else:
            integrated_error = error_over_time[0] * time_step  # Simple case for single time point

        print(f"Sum of squared differences: {sum_of_squares}")
        print(f"Integrated error over time: {integrated_error}")

        # Plot the error over time
        plt.figure(figsize=(10, 6))
        # plt.plot(time_points, error_over_time, 'r-', linewidth=2)
        plt.xlabel('Time')
        plt.ylabel('Squared Error')
        plt.title('Error Evolution Over Time')
        plt.grid(True)
        plt.savefig('error_over_time.png')
        plt.close()

        print("Plot saved as 'error_over_time.png'")

        # return sum_of_squares, integrated_error

    except Exception as e:
        print(f"Error: {e}")
        return None, None

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py file1.csv file2.csv [time_step]")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    # Optional time step parameter
    time_step = float(sys.argv[3]) if len(sys.argv) >= 4 else 1.0

    analyze_power_loss_error(file1_path, file2_path, time_step)
