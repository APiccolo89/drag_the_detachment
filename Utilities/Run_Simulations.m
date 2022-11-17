%========================================================================%
% Function utilities                                                     %
%========================================================================%
function [Tests] = Run_Simulations(D0_v,l0_v,eta0DS_v,DfS_v,n,s0_v,eta0DM_v,DfM_v,Benchmark,nlm)
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
if nlm.islinear 
    DfM_v = 0; 
end
[ium,iut0,iuet0,iuL0,iuDfS,iuDfM,iuD0]=ndgrid(eta0DM_v,s0_v,eta0DS_v,l0_v,DfS_v,DfM_v,D0_v);
tic
n_tests = length(ium(:));
disp(['Number of tests is ',num2str(n_tests)])
WORK = 1; 
for i=1:n_tests
    T_name   = strcat('T_',num2str(i));
    eta0DM   = ium(i);
    s0       = iut0(i);
    eta0DS   = iuet0(i);
    l0       = iuL0(i);
    Df_S     = iuDfS(i);
    Df_UM    = iuDfM(i);
    D0       = iuD0(i);
    [Temp]   = Processing_simulation(eta0DS,Df_S,n,l0,s0,D0,eta0DM,Benchmark,Df_UM,nlm);
    Tests.(T_name)=Temp;
    percentages_test = (i/n_tests)*100;
    disp([num2str(percentages_test,2),'per cent completed'])
end
toc
bla=1.0;
end


function [Temp]=Processing_simulation(eta0DS,Df_S,n,l0,s0,D0,eta0DM,Benchmark,Df_UM,nlm)

    [ID] = Compute_slab_characteristics(eta0DS,Df_S,n,l0,s0,D0,eta0DM,Df_UM,nlm);
    % Run a simulation with a specific combination of parameter
    Testdata = Run_Simulation_Drag(ID,Benchmark,nlm);
    Testdata_a = Run_Simulation_DragA(ID.ID_A,nlm);
    %interpolation routine from AD->D
    Temp = Testdata;
    Temp.Testdata_a=Testdata_a;
    Temp.initial_data = ID;
    ID = [];
    Testdata = [];
end

