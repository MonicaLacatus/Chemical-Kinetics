# Chemical Kinetics
 This repository contains several python files related to the field of reaction kinetics modelling.
### Reaction Equilibirum Assumption (REA):
The following code aims to demonstrate the use of the reaction equilibirum assumption (REA) by applying it
to a series reaction taking place inside a batch reactor: (A 🡘 B 🡘 C).
This is an assumption often used to simpliy reaction models and reduce the number of ODE balances needed to be solved. In the case of the above series reaction if k<sub>2</sub>, k<sub>-2</sub> >> k<sub>1</sub>, k<sub>-1</sub>  meaning that 
the second reaction is fast enough to equilibrate  immediately after any displacement from its equilibirum condition, 
than the full reaction model can be simplified on the basis that rate of the second reaction B 🡘 C is equal to zero. A comparison of the REA approximation with the full reaction model is supplied in this code.

### Quasi-Steady-State Assumption (QSSA):
The following code aims to demonstrate the use of the quasi-steady-state assumption (QSSA) by applying it
to a series reaction taking place inside a batch reactor: (A → B → C).
This is an assumption often used to simpliy reaction models. In the case of the series reaction in this example, if B is a highly reactive intermediate which is able to equilibrate quickly to its quasi-steady state value, its production rate can be set equal to zero. This simplifies the chemical reaction model significantly.A comparison of the QSSA approximation with the full reaction model is supplied in this code.
