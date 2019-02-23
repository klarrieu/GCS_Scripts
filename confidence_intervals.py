import numpy as np
import scipy.stats as stats


def r_to_z(r):
    return 0.5 * np.log((1+r)*1.0/(1-r))


def z_to_r(z):
    return (np.exp(2*z)-1)*1.0/(np.exp(2*z) + 1)


def r_confidence_interval(r, n, alpha=0.05):
    # returns confidence interval at the 1-alpha level for correlation of r with n observations
    # when alpha=0.05, it returns the range of possible population correlations at the 95% confidence level
    # so if 0 is not within the bounds, then the correlation is statistically significant at the 95% level
    z = r_to_z(r)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)  # 2-tailed z critical value

    lo = z - z_crit * se
    hi = z + z_crit * se

    # Return a sequence
    return (z_to_r(lo), z_to_r(hi))

print(r_confidence_interval(0.07, 1000, alpha=0.05))
print(r_confidence_interval(0, 1000, alpha=0.05))

def white_noise_confidence_interval(n):
    return (-1.0/n - 2.0/np.sqrt(n), -1.0/n + 2.0/np.sqrt(n))

def ar1_confidence_interval():
