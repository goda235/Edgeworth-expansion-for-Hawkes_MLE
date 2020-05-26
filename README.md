# Edgeworth expansion for Hawkes_MLE
By using the Monte Carlo method, numerically calculate the density function for the error of the maximum likelihood estimator of 1-dim exponential Hawkes process by 2nd order Edgeworth expansion.

Please execute the program in the order of the following numbers.

1. Hawkes_simulation.R
2. Hawkes_moments.R
3. Edgeworth_exp.R
4. Hawkes_asymptotic_normal.R
5. Q-Qplot.R

## 1.Hawkes_simulation.R
Simulate 1-dim Hawkes process by Ogata's method. Moreover, plot a path.

## 2. Hawkes_moments.R
Compute the expectation value of the derivative of the log-likelihood process numerically by the Monte Carlo method. 

## 3. Edgeworth_exp.R
Compute the density function of Edgeworth expansion numerically.

## 4. Hawkes_asymptotic_normal.R
The simulation is performed and a histogram is drawn.

## 5. Q-Qplot.R
Draw a Q-Q plot by using the simulation data from "4. Hawkes_asymptotic_normal.R".
