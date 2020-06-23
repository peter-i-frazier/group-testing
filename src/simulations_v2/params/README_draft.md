## Parameters for Stochastic Population-level Simulation

This folder contains files of parameter configurations under different scenarios. 

The file `fall_nominal.yaml` specifies the nominal parameter values for the following aspects of the model:
* population size
* initial prevalence
* duration of different stages of the disease
* asymptomatic fraction
* contact rates and disease transmission rates
* outside infection rate
* self-reporting rate
* contact tracing
* asymptomatic surveillance

Our model divides the population into five age groups. The file `age_severity_config.yaml` specifies the parameter values for
* the age distribution on Cornell's campus in the fall semester
* the probability of infection upon a close contact for different age groups
* the distribution of four levels of severity for an infected person in each age group. The four levels correspond to "asymptomatic", "mildly symptomatic but does not require hospitalization", "requires hospitalization but not critical care", and "requires critical care".

Other files in this folder import the `fall_nominal` file and makes modifications depending on the specific scenario.

(explain the june_8 folder and multigroup folder?)
