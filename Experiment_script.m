%==========================================================================
clear all
close all
clf
%==========================================================================
addpath Adimensional\
addpath Dimensional\
addpath Utilities\
addpath ../../../export_fig/
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
% Some important note: The numerical code has not a default viscosity cut
% off to avoid the limitation that usual 2D numerical code has. On the
% other hand, there are some combination of parameter that are not working
% well together, and this is a consequence that the non linear viscosity of
% of the slab is too unconstrained. xiUS must be used carefully, due to
% these reason.
%==========================================================================
% Folder output:
ptsave = ['../Tests_Results_LINEAR'];
% if the folder does not exist, create the folder.
if not(isdir(ptsave))
    mkdir(ptsave)
end
nlm = Problem_type;
nlm.Linear=0;   % Switching the position of linear-non_linear activate the non linear upper mantle routine.
nlm.iteration = 1;
nlm.cut_off   = 0;

if nlm.islinear == 1
    ptsave = ['../Tests_Results_LINEAR'];

else
    ptsave = ['../Tests_Results_NLINEAR'];

end 
xiUM_v = [];
Limit_scenario = 'Linear_Regime_Main';

switch Limit_scenario
    case 'Linear_Regime_Main' % xim < 1e-6
        %Geometric properties of the Slab
        l0_v      = (300e3:50e3:600e3);         % initial length [m]
        D0_v      = 80e3;                      % thickness [m]
        s0_v      = [60e6:30e6:240e6];        % reference buoyancy stress [Pa]
        % Slab Rheology
        eta0DS_v   = 10.^[22,24,26,28,30];                      % [Pas] reference diffusion creep viscosity of the slab
        xiUS_v     = 10.^[-1.0,0.0,1.0,2.0,6.0];              % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress
        n_v        = [2.0,3.5,7.0];            % [n.d.] pre-exponential factor
        % Upper Mantle
        eta0DM_v  = 10.^[19,20,21,22,23];                       % [Pas] reference diffusion creep viscosity of the upper mantle
        xiUM_v    = 10.^[-12,-10,-8,-6,-4];                  % [n.d.] viscosity contrast between diffusion and dislocation creep at reference stress (UM)

    case 'Transition_Regime_Main' % xim >1e-6 & xim<1e6
        l0_v      = (300e3:50e3:600e3);         % initial length [m]
        D0_v      = 80e3;                      % thickness [m]
        s0_v      = [60e6:30e6:240e6];        % reference buoyancy stress [Pa]
        % Slab Rheology
        eta0DS_v   = 10.^[22,24,26,28,30];                      % [Pas] reference diffusion creep viscosity of the slab
        xiUS_v     = 10.^[-1.0,0.0,1.0,2.0,6.0];              % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress
        n_v        = [3.5];            % [n.d.] pre-exponential factor
        % Upper Mantle
        eta0DM_v  = 10.^[19,20,21,22,23];                       % [Pas] reference diffusion creep viscosity of the upper mantle
        xiUM_v    = 10.^[-2,-1,0,1,2];                  % [n.d.] viscosity contrast between diffusion and dislocation creep at reference stress (UM)

    case 'NonLinear_Regime_Main'
        l0_v      = (300e3:50e3:600e3);         % initial length [m]
        D0_v      = 80e3;                      % thickness [m]
        s0_v      = [60e6:30e6:240e6];        % reference buoyancy stress [Pa]
        % Slab Rheology
        eta0DS_v   = 10.^[22,24,26,28,30];                      % [Pas] reference diffusion creep viscosity of the slab
        xiUS_v     = 10.^[-1.0,0.0,1.0,2.0,6.0];              % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress
        n_v        = [3.5];            % [n.d.] pre-exponential factor
        % Upper Mantle
        eta0DM_v  = 10.^[19,20,21,22,23];                       % [Pas] reference diffusion creep viscosity of the upper mantle
        xiUM_v    = 10.^[4,6,8,10,12];                  % [n.d.] viscosity contrast between diffusion and dislocation creep at reference stress (UM)
end
%==========================================================================

for i=1:length(n_v)
    n = n_v(i);
    name_tests= strcat('Tests_n_',num2str(n));
    % Function to run the ensamble of test
    [Tests] = Run_Simulations(D0_v,l0_v,eta0DS_v,xiUS_v,n,s0_v,eta0DM_v,xiUM_v,Benchmark,nlm);
    % Update the test structure
    T.(strcat('n_',num2str(floor(n))))= Tests;
    Tests = [];
end
%% Save information into DataBase
[Data_S] = extract_information_detachment(T,1,nlm);

name_data_base = strcat('../Data_Base/',Limit_scenario,'_data_base_gamma_1000e3_exp','.mat');

initial_vectors.D0_v = D0_v;
initial_vectors.eta0DM_v = eta0DM_v;
initial_vectors.eta0DS_v=eta0DS_v;
initial_vectors.l0_v  = l0_v;
initial_vectors.n_v = n_v;
initial_vectors.xiUM_v = xiUM_v;
initial_vectors.xiUS_v = xiUS_v;
initial_vectors.s0_v   = s0_v;
if nlm.islinear == 1
    filename_ = '../Data_Base/Linear_Tests_Data_Base.mat';
    xiUM_v = 0;
    [nTv,~] = ndgrid(eta0DM_v,s0_v,eta0DS_v,l0_v,xiUS_v,xiUM_v,D0_v,n_v);
    number_tests = length(nTv(:));
    initial_vectors.number_tests = number_tests;
else
    filename_ = name_data_base;
    [nTv,~] = ndgrid(eta0DM_v,s0_v,eta0DS_v,l0_v,xiUS_v,xiUM_v,D0_v,n_v);
    number_tests = length(nTv(:));
    initial_vectors.number_tests = number_tests;
end
n_3 = T.n_3;
save(filename_,'Data_S','n_3','initial_vectors')
