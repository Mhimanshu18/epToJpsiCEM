import pandas as pd
import matplotlib.pyplot as plt

# Function to read the data from a .dat file and return the x and y columns
def read_data(file_name):
    data = pd.read_csv(file_name, delim_whitespace=True, header=None)
    x = data[0]  # First column (x values)
    y = data[1]  # Second column (y values)
    return x, y

# File names
files = [
    'data/gmmr_2012_158.1.dat', 'data/gmmr_2012_31.6.dat', 'data/icem_158.1.dat', 'data/icem_31.6.dat',
    'data/gmmr_2012_17.33.dat', 'data/gmmr_2012_7.2.dat', 'data/icem_17.33.dat', 'data/icem_7.2.dat'
]

# Set up the plot
plt.figure(figsize=(10, 6))

# Loop through each file, read the data, and plot it
for file in files:
    x, y = read_data(file)
    label = file.split('.')[0]  # Use the file name (without extension) as the label
    plt.plot(x, y, label=label)

# Customize the plot
plt.title('Comparison of Data from Different Files')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.legend(loc='best')
plt.grid(True)

# Show the plot
plt.show()
