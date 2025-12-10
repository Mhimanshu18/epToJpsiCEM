import numpy as np
import matplotlib.pyplot as plt

# Load data from the .dat file
# Update the filename and delimiter based on your file format
#filename = 'gmmr_2012_7.2.dat'
data = np.loadtxt('data/gmmr_2012_7.2.dat', delimiter='\t')  # Change delimiter if needed (e.g., '\t' for tab-separated)

# Assuming two columns: X and Y
x = data[:, 0]  # First column
y = data[:, 1]  # Second column

# Plotting the data
plt.figure(figsize=(8, 6))
plt.plot(x, y, label='CEM', marker='o', linestyle='-', color='b')

# Adding labels, title, and legend
plt.xlabel('qt', fontsize=12)
plt.ylabel('Asymmetry', fontsize=12)
plt.title('Asymmetry vs qt', fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout()

# Show or save the plot
plt.savefig('epToJpsiCEM.png')  # To save the plot
plt.show()
