This page includes the R and BUGS code for the tutorial paper on joint modeling as described in Baghfalaki et al. (2024). [https://arxiv.org/abs/2403.07778]

#### Summary
An introduction to the principles of Bayesian joint modeling of longitudinal measurements and time-to-event outcomes, as well as model implementation using the BUGS language syntax. This syntax can be executed directly using OpenBUGS or by using convenient functions to call OpenBUGS and JAGS from R software. All the details of joint models are provided, ranging from simple to more advanced models. The discussion began with the joint modeling of a Gaussian longitudinal marker and time-to-event outcome. The implementation of the Bayesian paradigm of the model is reviewed. The strategies for simulating data from the joint model are also discussed. The details of the paper are provided in Baghfalaki et al. (2024), Joint modeling of longitudinal measurements and time-to-event outcomes using BUGS. 

This page includes the R and BUGS code for the following joint modeling. For more details about this part, refer to Baghfalaki et al. (2024).
### JM FOR ONE GAUSSIAN LONGITUDINAL MARKER
#### Simulating data for JM with a Gaussian longitudinal marker 
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/Simulated_data.

#### JM with a proportional hazard sub-model with constant baseline hazard
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/constant_BH.

#### JM with a proportional hazard sub-model with a Weibull baseline hazard
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/Weibull_BH.

#### JM with a proportional hazard sub-model with a piecewise constant baseline hazard
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/PC_BH.

#### JM with a proportional hazard sub-model with a spline baseline hazard
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/BS_PH.

#### Considering current value and current slope as association between two sub-models
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/current_slope.

#### A shared random effects model for the association between the two sub-models
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/BS_PH.

### JM FOR MULTIVARIATE LONGITUDINAL MARKERS
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/multivariate.

### JM FOR LONGITUDINAL MEASUREMENTS AND COMPETING RISKS OUTCOME
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/CR.

### JM OF ZERO-INFLATED LONGITUDINAL MARKER AND SURVIVAL TIME
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/Zero_inflation.


### JM OF LONGITUDINAL MARKER AND SURVIVAL TIME with a CURE FRACTION 
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/Applications/JMC

### REAL DATA ANALYSIS
 The codes can be found at https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/Applications


### Reference
Baghfalaki, T., Ganjali, M., Barbieri, A., Hashemi, R. and Jacqmin-Gadda, H. (2024). Joint Modeling of Longitudinal Measurements and Time-to-event Outcomes Using BUGS, https://arxiv.org/abs/2403.07778.


