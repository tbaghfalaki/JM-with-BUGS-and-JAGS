
# Simulated data
Here, we will examine two distinct simulated datasets. 

## 1
The first one is based on a constant value for the baseline hazard, so the inverse of the cumulative hazard function is a common method for this purpose. The R code can be found in "simulated_JM1.R", and the simulated data can be found in "long.data_1.RData" and "surv.data_1.RData".
The longitudinal profiles and KM plot are presented [here](/Figures/sim.md). 

## 2 Example 2: Generation using the permutational algorithm
The real values of the parameters are the same as the first set, but here a uniform distribution is considered for generating censoring time. The data generation is based on the permutational algorithm using the PermAlgo package. The R code for this method can be found in "simulated_JM2.R", and the simulated data can be found in "long.data_2.RData" and "surv.data_2.RData". 
The longitudinal profiles and KM plot are presented [here](/Figures/sim.md).
