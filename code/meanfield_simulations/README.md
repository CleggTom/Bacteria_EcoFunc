# Meanfield Approaches

This folder contains the julia and R code to simulate and analyse the meanfield approximations. The files are as follows:

##Random Matrix
- `temp_simulations.jl`: Simulation of a communities across a temperature gradient. simulations are replicated for each temperature with randomly assembled communities.

- `temp_meanfield_analysis.R`: Analysis of the temperature simulations showing the quality of the approximations.

- `interaction_simulations.jl`: Does the same as `temp_simulations.jl` but varies the average interaction strength across a gradient

- `int_meanfield_analysis.R` Analysis of the varying interaction strength parameters.

##Symmetric Matrix
