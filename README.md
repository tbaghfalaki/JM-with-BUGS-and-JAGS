# Joint modeling of longitudinal measurements and time-to-event outcomes using BUGS
### Taban Baghfalaki | Mojtaba Ganjali | Antoine Barbieri | Reza Hashemi | Hélène Jacqmin-Gadda
An introduction to the principles of Bayesian joint modeling of longitudinal measurements and time-to-event outcomes, as well as model implementation using the BUGS language syntax. This syntax can be executed directly using OpenBUGS or by using convenient functions to call OpenBUGS and JAGS from R software. All the details of joint models are provided, ranging from simple to more advanced models. The discussion began with the joint modeling of a Gaussian longitudinal marker and time-to-event outcome. The implementation of the Bayesian paradigm of the model is reviewed. The strategies for simulating data from the joint model are also discussed. The details of the paper are provided in Baghfalaki et al. (2024), Joint modeling of longitudinal measurements and time-to-event outcomes using BUGS. 

This page includes the R and BUGS code for the following joint modeling. 

### JM FOR ONE GAUSSIAN LONGITUDINAL MARKER
#### Simulating data for JM with a Gaussian longitudinal marker 
For this (https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/BS_PH)

#### JM with a proportional hazard sub-model with constant baseline hazard

#### JM with a proportional hazard sub-model with a Weibull baseline hazard

#### JM with a proportional hazard sub-model with a piecewise constant baseline hazard

#### JM with a proportional hazard sub-model with a spline baseline hazard

#### Considering current value and current slope as association between two sub-models


#### A shared random effects model for the association between the two sub-models

### JM FOR MULTIVARIATE LONGITUDINAL MARKERS

### JM FOR LONGITUDINAL MEASUREMENTS AND COMPETING RISKS OUTCOME

### JM OF ZERO-INFLATED LONGITUDINAL MARKER AND SURVIVAL TIME

### REAL DATA ANALYSIS


####
