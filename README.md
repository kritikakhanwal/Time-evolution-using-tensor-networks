# Title : Time-evolution-using-tensor-networks
Time evolution using TEBD method. 

# Installation: 
This project extensively uses the amazing library in Julia [ITensor](https://itensor.github.io/ITensors.jl/dev/). 

To proceed first Julia installation is needed and then install ITensor library using import Pkg; Pkg.add("ITensors")

# Description:
The purpose of this project is to get some hands-on experience with Tensor structures that can be used for simulations of quantum systems. Further, we have time evolved the product state and also found the entanglement entropy of the state. The time evolution is done using TEBD method. 

We also tried to time evolve a density matrix starting from an identity (maximally entangled state) using the same TEBD method. This is interesting since the density matrix will be a matrix-product operator. However, we can write it as an MPS and then do the time evolution of the purified state. The method is explained in pdf. More details are presented in this nice [paper](https://arxiv.org/abs/1910.09142). 

