# Edgeworth expansion for Hawkes_MLE
By using the Monte Carlo method, numerically calculate the density function for the error of the maximum likelihood estimator of the 1-dim exponential Hawkes process by 2nd order Edgeworth expansion.
Finally, draw a histgram and Q-Q plot.

## Setting
params  : Hawkes parameter

t_max : Observed time

MC_coef : Number of the Monte Carlo simulation for computing coefficients of Edgeworth expansion

MC_boot : Number of the Monte Carlo simulation for computing Bootstrap density

MC_hist : Number of the Monte Carlo simulation for computing MLEs

## Main Functions
### simulate_uni_hawkes(params, t_max)
Simulate 1-dim Hawkes process by Ogata's method.

### loglik(params, arrivals, t_max)
The log-likelihood function of the 1-dim exponential Hawkes process. "arrivals" is a data of Hawkes process.

### Coefficients(params, t_max, MC_coef)
Compute the density function of Edgeworth expansion numerically.

### MC_MLE(params, MC)
Compute MC number of MLEs.
