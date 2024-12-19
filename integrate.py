import scipy.integrate as integrate

# Define the function y(x) that you want to integrate
def y(x):
    return x  # Example function: y(x) = x^2

# Perform the integration over the range from 0 to 1
integral, error = integrate.quad(y, 0, 1)

# Print the result of the integration
print(f"Integral of y(x) from 0 to 1: {integral}")
print(f"Estimated error in the integration: {error}")
