# mig-kernel

Code for the paper 

> Normal approximations for the multivariate inverse Gaussian distribution and
asymmetric kernel smoothing on $d$-dimensional half-spaces

by Alain Desgagné, Christian Genest, Frédéric Ouimet and Léo R. Belzile

The corresponding **R** package can be found on Github at: https://github.com/lbelzile/mig

**Files**:

- `simulation_study_MIG_R1.R` launches the simulation study. Scenarios were run sequentially using the `simsalapar` package in **R** on a local computer.
- `utils.R` contains utilities to generate from the simulation scenarios and evaluate true density for the four distributions $$F_1--F_4$$. It contains the workhorse `doOne` that generates results in the simulation study for a single scenario.
-
