import numpy as np
import matplotlib.pyplot as plt

# Load data from the .dat files
cem = np.loadtxt('data/gmmr_qT_2012_158.1.dat', delimiter='\t')  # Update path and delimiter as needed
icem = np.loadtxt('data/icem_qT1_158.1.dat', delimiter='\t')

# Extract X (first column) and Y (second column) for both datasets
cem_x = cem[:, 0]
cem_y = cem[:, 1]
icem_x = icem[:, 0]
icem_y = icem[:, 1]

# Create a subplot for comparison
fig, axs = plt.subplots(1, 1, figsize=(10, 6))  # Single plot for comparison

# Plot cem and icem data
axs.plot(cem_x, cem_y, label='cem_158.1', color='blue', marker='o', linestyle='-')
axs.plot(icem_x, icem_y, label='icem_158.1', color='red', marker='s', linestyle='--')
# Customize the plot
axs.set_title('Comparison: cem_158.1 vs icem_158.1', fontsize=14)
axs.set_xlabel('qT', fontsize=12)
axs.set_ylabel('Asymmetry', fontsize=12)
axs.legend(fontsize=10)
axs.grid(True)

# Adjust layout and save/show the plot
plt.tight_layout()
plt.savefig('comparison_qT_158.1.png')  # Save the plot as an image
plt.show()
