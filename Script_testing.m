%==========================================================================
clear all
close all
clf
%==========================================================================
addpath Adimensional\
addpath Dimensional\
addpath Utilities\
%==========================================================================
%Parameter to test: 
% To do: 
% Generate a data structure and labelling. 

l0_v    = (300e3:50e3:600e3);                 % initial length
D0      = 80e3;                  % thickness

% Slab Rheology
eta0_v = [1e22,5e22,1e23,5e23];                   % [Pas] refernce power law viscosity slab 
Df_v   = [0.8,1.0,10,100];                          % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress 
n_v     =[2.0, 3.5,7.0,9.0];                       % [n.d.] pre-exponential factor
s0_v  = [100e6:50e6:200e6];                   % [Pa]  reference buoyancy stress

% Upper Mantle
etaum_v =10.^(18:0.1:21);    % [Pa.s]vector of the mantle viscosity

name_tests_v = {'D0_80km_n_2dot0_Df_EXT','D0_80km_n_3dot5_Df_EXT','D0_80km_n_7dot0_Df_EXT','D0_80km_n_9dot0_Df_EXT'};
ptsave = 'Tests_result';
if not(isfolder(ptsave))
    mkdir(ptsave)
end

for i=1:4
    name_tests = name_tests_v{i};
    n = n_v(i);

    % Function to run the ensamble of test 
    [Tests] = Run_Simulations(D0,l0_v,eta0_v,Df_v,n,s0_v,etaum_v);

    plot_results(Tests,name_tests,ptsave)
end
% Function to plot the results 
