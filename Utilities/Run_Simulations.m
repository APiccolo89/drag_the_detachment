%========================================================================%
% Function utilities                                                     %
%========================================================================%
function [Tests] = Run_Simulations(D0,l0_v,eta0_v,Df_v,n,s0_v,etaum_v)
%=====================================================================
% Output:
% ====================================================================
% [Testdata] => structure containing the results of each of the
%               simulation performed
%=====================================================================
% Input:
%
% D0,l0_v => Initial thickness and length (l0_v vector of data)
%
% eta0_v => reference power law viscosity at reference stress (tau0)
%       n power law exponent [vector]
%
% Df    =>viscosity contrast between diffusion/dislocation creep @
%       reference stress (tau0)
%
% etaum_v => Viscosity of the upper mantle [vector].
%=====================================================================
% Alg:
% -> loop over the data and generate combination of simulation
% -> construct the characteristic values
% -> run simulation per each combination of parameter
%=====================================================================
% create a multidimensional array containing all the initial data
[ium,iut0,iuet0,iuL0,iuDf]=ndgrid(etaum_v,s0_v,eta0_v,l0_v,Df_v);
tic
n_tests = length(ium(:));
for i=1:n_tests
    etaum = ium(i);
    s0    = iut0(i);
    eta0  = iuet0(i);
    l0    = iuL0(i);
    Df    = iuDf(i);
    [ID] = Compute_slab_characteristics(eta0,Df,n,l0,s0,D0,etaum);
    % Create a labeling for the simulation
    str0   = strcat('Sim_D0_',num2str(int16(D0/1e3)),num2str(i));
    l_simulation = strcat(str0);
    disp(['Simulation ',l_simulation, 'is starting'])
    % Run a simulation with a specific combination of parameter
    Testdata = Run_Simulation_Drag(ID);
    %Testdata_a = Run_Simulation_DragA(ID.ID_A);
    Tests.(l_simulation) = Testdata;
    %Tests.(l_simulation).Testdata.A=Testdata_a;
    Tests.(l_simulation).initial_data = ID;
    ID = [];
    Testdata = [];
    percentages_test = (i/n_tests)*100;
    disp([num2str(percentages_test,2),'per cent completed'])
end
toc
end

