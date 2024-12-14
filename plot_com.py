import numpy as np
import matplotlib.pyplot as plt

# Load the data for the four files
icem_158_1 = np.loadtxt("data/icem_158.1.dat")
gmmr_158_1 = np.loadtxt("data/gmmr_2012_158.1.dat")
icem_31_6 = np.loadtxt("data/icem_31.6.dat")
gmmr_31_6 = np.loadtxt("data/gmmr_2012_31.6.dat")
icem_17_33 = np.loadtxt("data/icem_17.33.dat")
gmmr_17_33 = np.loadtxt("data/gmmr_2012_17.33.dat")
icem_7_2 = np.loadtxt("data/icem_7.2.dat")
gmmr_7_2 = np.loadtxt("data/gmmr_2012_7.2.dat")

# Create a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# First comparison: icem_158.1 vs gmmr_158.1
axs[0, 0].plot(icem_158_1[:, 0], icem_158_1[:, 1], label='icem_158.1', color='blue')
axs[0, 0].plot(gmmr_158_1[:, 0], gmmr_158_1[:, 1], label='gmmr_158.1', color='red')
axs[0, 0].set_title('Comparison: icem_158.1 vs gmmr_158.1')
axs[0, 0].set_xlabel('Y')
axs[0, 0].set_ylabel('Asymmetry')
axs[0, 0].legend()

# Second comparison: icem_31.6 vs gmmr_31.6
axs[0, 1].plot(icem_31_6[:, 0], icem_31_6[:, 1], label='icem_31.6', color='blue')
axs[0, 1].plot(gmmr_31_6[:, 0], gmmr_31_6[:, 1], label='gmmr_31.6', color='red')
axs[0, 1].set_title('Comparison: icem_31.6 vs gmmr_31.6')
axs[0, 1].set_xlabel('Y')
axs[0, 1].set_ylabel('Asymmetry')
axs[0, 1].legend()

# Third comparison: icem_17.33 vs gmmr_17.33
axs[1, 0].plot(icem_17_33[:, 0], icem_17_33[:, 1], label='icem_17.33', color='blue')
axs[1, 0].plot(gmmr_17_33[:, 0], gmmr_17_33[:, 1], label='gmmr_17.33', color='red')
axs[1, 0].set_title('Comparison: icem_17.33 vs gmmr_17.33')
axs[1, 0].set_xlabel('Y')
axs[1, 0].set_ylabel('Asymmetry')
axs[1, 0].legend()

# Fourth comparison: icem_7.2 vs gmmr_7.2
axs[1, 1].plot(icem_7_2[:, 0], icem_7_2[:, 1], label='icem_7.2', color='blue')
axs[1, 1].plot(gmmr_7_2[:, 0], gmmr_7_2[:, 1], label='gmmr_7.2', color='red')
axs[1, 1].set_title('Comparison: icem_7.2 vs gmmr_7.2')
axs[1, 1].set_xlabel('Y')
axs[1, 1].set_ylabel('Asymmetry')
axs[1, 1].legend()

# Adjust the layout to prevent overlap
plt.tight_layout()

# Save the plot to a PDF file
plt.savefig("comparison_plots.pdf", format="pdf", dpi=300)

# Show the plot
plt.show()
