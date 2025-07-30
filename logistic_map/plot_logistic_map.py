import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = False  # Ensure that external LaTeX is not used
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'  # Even if LaTeX is not used, this can help with math rendering
plt.rcParams['mathtext.fontset'] = 'cm'  # Use the Computer Modern font for math text

# Define the range of values for r
r_values = np.arange(3, 4.0, 0.0001)

# Number of iterations and the number to plot
iterations = 2000
last = 1000

# Initialize the array for x values
x = 1e-5 * np.ones(r_values.shape)

# Create the figure and axis
plt.figure(figsize=(10, 6))

# Iterate over r values
for i in range(iterations):
    x = r_values * x * (1 - x)
    # For the last iterations, plot the points
    if i >= (iterations - last):
        plt.plot(r_values, x, ',k', alpha=0.75, markersize=0.5)

# Configure the axes
plt.xlim(3, 4)
plt.title("Bifurcation diagram of the logistic map")
plt.xlabel(r"$\mu$")
plt.ylabel("x")
plt.show()
