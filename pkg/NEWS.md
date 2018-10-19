# pulseR 1.0.3

* `initParams` assigns the first normalisation factor to 1 if
  the normalisation factors are used
* update vignettes, add a vignette for the RNA modification to model m6A

# pulseR 1.0.2

* update documentation
* the stopping criteria is based on the **ABSOLUTE** change of 
parameters, not relative. If the the parameter can have zero value,
relative changes make no sense. It is advised to provide parameters
on the same scale then, e.g. log(expression level) instead of
working the mean read counts units.

# pulseR 1.0.1

* fixed order of formulae in PulseData object. Now they are sorted according
  to the factor levels in the condition data.frame
  
# pulseR 1.0.0

* Added confidence interval estimations



