%==========================================================================
clear all
close all
clf
%==========================================================================
addpath Adimensional\
addpath Dimensional\
addpath Utilities\
Benchmark = 0.0 ; % Benchmark activaction flag
%==========================================================================
%Parameter to test: 
%==========================================================================
% Testing parameter: 
%==========================================================================
% l0_v = vector containing all the lenght to test
% D0   = initial thickness 
% s0_v = buoyancy stress applied 
%==========================================================================
% eta0DS_v = [vector] reference linear viscosity of the slab.
% eta0DM_v = [vector] reference linear viscosity of the mantle.
% DfS_v    = [vector] viscosity contrast between power-law rheology and
% linear one at reference stress condition
% DfM_v    = [vector] viscosity contrast between power-law rheology and
% linear one at reference condition [deactivated if upper mantle is linear]
% n_v      = [vector] power law exponent {same Slab-Mantle}
% nlm      = activation flag non linear mantle. 
% These are the primary values used to construct the set of experiments.
% Instead of using directly Lambda as primary value, I decided to use
% Earth-like range parameters that are used to compute the relevant data
% such that is possible to associated set of parameters to the
% dimensionless group. 
%==========================================================================
%[OUTPUT]: => [Tests.mat] Database featuring relevant data for each tests
%i.e. Table_Experiments:     T1 => [D1-------Dn|Detachment data]
%     Numerical_Experiments: T1 -
%                               |_ INITIAL DATA STRUCTURE
%                               |_ D_norm [1,nstep] 
%                               |_ Stress [3,nstep]: 1:Effective stress 2:
%                                   Bouyancy stress 3: Drag integral stress
%                               |_ Def [3, nstep]:  1: total strain; 2:
%                                   dislocation creep strain; 3: Diffusion
%                                   creep strain
%                               |_ Lambda[1,nstep]
%                               |_ UM_vis[1,nstep]
%                               |_ time_Vector[1,nstep]
%                            Tn    
% The numbering of the test depends on the layout of the matrix created in 
% f. Run Simulations[]
% Pictures: saved outside the actual folder of the script. And divided
% using the stress exponent.
%==========================================================================
% Folder output: 
ptsave = '../Tests_Results'; 
% if the folder does not exist, create the folder. 
if not(isfolder(ptsave))
    mkdir(ptsave)
end

%Geometric properties of the Slab 
l0_v      = (300e3:10e3:600e3);         % initial length [m]
D0_v      = 80e3;                       % thickness [m]
s0_v      = [100e6:100e6:200e6];        % reference buoyancy stress [Pa]
% Slab Rheology
eta0DS_v  = [1e22,1e23,1e24];                      % [Pas] reference diffusion creep viscosity of the slab 
DfS_v     = [1.0, 10, 100];               % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress 
n_v       = [ 3.5,2.0,7.0,9.0];            % [n.d.] pre-exponential factor
% Upper Mantle
eta0DM_v  = [1e19,1e20,1e21,1e22];                       % [Pas] reference diffusion creep viscosity of the upper mantle 
DfM_v     = 10.^(3:7);                  % [n.d.] viscosity contrast between diffusion and dislocation creep at reference stress (UM)
nlm       = Problem_type.NonLinear;     % Switching the position of linear-non_linear activate the non linear upper mantle routine. 
%==========================================================================

for i=1:4
    n = n_v(i);
    name_tests= strcat('Tests_n_',num2str(n));
    % Function to run the ensamble of test 
    [Tests] = Run_Simulations(D0_v,l0_v,eta0DS_v,DfS_v,n,s0_v,eta0DM_v,DfM_v,Benchmark,nlm);
    plot_results(Tests,name_tests,ptsave)
end

