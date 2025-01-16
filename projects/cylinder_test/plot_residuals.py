import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

def plot_newton_residuals(filename):
    # Load data, skip header
    data = np.loadtxt(filename, skiprows=1)
    timesteps = data[:, 0]
    iterations = data[:, 1]
    residuals = data[:, 2]
    times = data[:, 3]
    currents = data[:, 4]

    # Create figure with three subplots
    fig = plt.figure(figsize=(15, 12))
    gs = GridSpec(3, 1, height_ratios=[2, 1, 1])
    
    # Plot 1: Residuals vs iterations for each timestep with colormap
    ax1 = fig.add_subplot(gs[0])
    unique_timesteps = np.unique(timesteps)
    
    # Create colormap
    cmap = plt.cm.viridis
    norm = Normalize(vmin=times.min(), vmax=times.max())
    
    for timestep in unique_timesteps:
        mask = timesteps == timestep
        current_time = times[mask][0]
        ax1.semilogy(iterations[mask], residuals[mask], 'o-', 
                    color=cmap(norm(current_time)))
    
    # Add colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm, ax=ax1)
    cbar.set_label('Time [s]')
    
    ax1.grid(True)
    ax1.set_xlabel('Newton Iteration')
    ax1.set_ylabel('Residual (log scale)')
    ax1.set_title('Newton Method Convergence per Timestep')
    
    # Plot 2: Number of Newton iterations per timestep
    ax2 = fig.add_subplot(gs[1])
    iterations_per_timestep = [sum(timesteps == t) for t in unique_timesteps]
    unique_times = times[np.unique(timesteps, return_index=True)[1]]
    ax2.plot(unique_times, iterations_per_timestep, 'b-o')
    ax2.grid(True)
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Newton Iterations')
    ax2.set_title('Newton Iterations Required per Timestep')
    
    # # Plot 3: Input current over time with convergence info
    # ax3 = fig.add_subplot(gs[2])
    # ax3.plot(unique_times, currents[np.unique(timesteps, return_index=True)[1]], 'r-', label='Current')
    
    # # Add markers for timesteps requiring many iterations
    # high_iter_mask = iterations_per_timestep > np.mean(iterations_per_timestep) + np.std(iterations_per_timestep)
    # if np.any(high_iter_mask):
    #     ax3.plot(unique_times[high_iter_mask], 
    #             currents[np.unique(timesteps, return_index=True)[1]][high_iter_mask],
    #             'k*', markersize=10, label='High iteration count')
    
    # ax3.grid(True)
    # ax3.set_xlabel('Time [s]')
    # ax3.set_ylabel('Current')
    # ax3.set_title('Input Current vs Time')
    # ax3.legend()
    
    plt.tight_layout()
    plt.savefig('newton_convergence_analysis.png', bbox_inches='tight', dpi=300)
    plt.close()

    # Print statistics
    print("\nConvergence Statistics:")
    print(f"Total number of timesteps: {len(unique_timesteps)}")
    print(f"Average iterations per timestep: {np.mean(iterations_per_timestep):.2f}")
    print(f"Standard deviation: {np.std(iterations_per_timestep):.2f}")
    print(f"Maximum iterations in a timestep: {np.max(iterations_per_timestep)}")
    print(f"Minimum iterations in a timestep: {np.min(iterations_per_timestep)}")
    print(f"\nWorst convergence at time = {unique_times[np.argmax(iterations_per_timestep)]:.3f}s")
    
if __name__ == "__main__":
    plot_newton_residuals("newton_residuals.txt")