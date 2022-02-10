# GBC_DOP_Letscher
Matlab code associated with Letscher et al 2022 in GBC

doit_DOP.m        -->  the main driver script

neglogpost_DOP.m  -->  computes the cost function

eqPcycle_DOP.m    -->  computes the equilibrium 3D DIP, DOP, POP distributions 
                       computes the first and second derivatives of the model solution w.r.t. the model parameters
                       
get_bottle_dom.m  -->  reads in the DOP observations from .csv file and interpolates to OCIM2 grid

buildPFD.m        -->  builds the sinking particle flux divergence operator

mfactor.m         -->  uses sparse LU to factorize a matrix or solve a linear system from prefactored matrices 

d0.m              -->  builds a sparse diagonal matrix given a 3D field


Other large data and operator files can be obtained from corresponding author robert.letscher@unh.edu
    e.g. OCIM2_CTL_He.mat, WOAPO4.mat, etc.
