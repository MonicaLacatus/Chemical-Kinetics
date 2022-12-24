# Chemical Kinetics
 This repository contains several python files related to the field of reaction kinetics modelling.
### Reaction Equilibirum Assumption (REA):
The following code aims to demonstrate the use of the reaction equilibirum assumption (REA) by applying it
to a series reaction taking place inside a batch reactor: (A 🡘 B 🡘 C).
This is an assumption often used to simpliy reaction models and reduce the number of ODE balances needed to be solved. In the case of the above series reaction if k<sub>2</sub>, k<sub>-2</sub> >> k<sub>1</sub>, k<sub>-1</sub>  meaning that 
the second reaction is fast enough to equilibrate  immediately after any displacement from its equilibirum condition, 
than the full reaction model can be simplified on the basis that rate of the second reaction B 🡘 C is equal to 0. A comparison of the REA approximation with the full reaction model is supplied in this code.
