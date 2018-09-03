# TransientGrowth
Project with UniMelb Aug 2018
Jimmy Phillip, Andre Peterson, Blake Cook

MATLAB code modified and maintained by Blake



Transient growth calculations for Neural Networks: Matrices of the form
given by Rajan and Abbott (PRL - 2006)

Transient growth in system given by 
  $$\frac{du_i}{dt} = -\frac{u_i}{tau} + \Sum_j W_{ij} \phi(u_j)$$



**** Requires MATLAB R2018a **** 

src/ - solvers and plotting functions
test/ - location of test files, re-run when major changes are made
output/ - location of all output

All code should be executed from the directory 'TransientGrowth', to access 
functions in 'src/', run 'addpath(genpath('src/'))' (similarly for 'test/')


Main file for demonstrating an individual run is 'TG_main.m'

