"""
Demonstration of astrophotographical stacking to
detect the underlying signal from an ensemble of
100s of samples.

"""
import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import norm

# Let's create a array of evenly spaced dense points covering 99% of the distribution.
x = np.linspace(norm.ppf(0.001), norm.ppf(0.999), 500, dtype=np.float32).reshape(-1,1)

# Now, let's find the pdf of the distribution at those points without adding any noise.
gaussian_without_noise = norm.pdf(x).reshape(-1,1)

# Let's plot our distribution
plt.plot(x, gaussian_without_noise)
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (m)")
plt.title("Pure Signal.")
plt.show()

# Now, we will create a list to hold 100 noisy distributions
noisy_dists = []

# Let's fill the list with noisy distributions
for i in range(100):
    noisy_dists.append(gaussian_without_noise + np.random.randn(*gaussian_without_noise.shape))

# Let's plot a couple of these noisy distributions to see how much
# noisy they are
fig, ax = plt.subplots(nrows=2, ncols=1)
ax[0].plot(x, noisy_dists[0])
ax[1].plot(x, noisy_dists[1])
fig.show()
plt.show()

# Now, let's mean stack them together to make a distribution with high SNR ratio
# so we can easily see the underlying signal.
noisy_dists = np.array(noisy_dists)
mean_stacked_dists = (1./100.)*np.sum(noisy_dists, axis=0)

# Let's plot our results
plt.plot(x, mean_stacked_dists, color = 'y', label="stacked_signal")
plt.plot(x, gaussian_without_noise, color = 'r', label="true signal")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (m)")
plt.title("Pure vs Stacked Signal.")
plt.legend()
plt.show()