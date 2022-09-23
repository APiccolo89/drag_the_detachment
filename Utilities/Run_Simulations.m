%========================================================================%
% Function utilities                                                     %
%========================================================================%
function [Tests] = Run_Simulations(D0,l0_v,eta0_v,Df_v,n,s0_v,etaum_v,Df_UM,Benchmark)
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
index_test=1:n_tests;
disp(['Number of tests is ',num2str(n_tests)])
% Create field structure
for i=1:n_tests
    str0   = strcat('Sim_D0_',num2str(int16(D0/1e3)),'_',num2str(index_test(i)));
    l_simulation = strcat(str0);
    disp(['Simulation ',l_simulation, 'is starting'])
    Tests.(l_simulation) = [];
    fields_structure{i}=l_simulation;
end



pause(5)
WORK = 1; 
for i=1:n_tests
    
    etaum = ium(i);
    s0    = iut0(i);
    eta0  = iuet0(i);
    l0    = iuL0(i);
    Df    = iuDf(i);
    [Temp]=Processing_simulation(eta0,Df,n,l0,s0,D0,etaum,Benchmark,Df_UM,fields_structure{i});
    Tests.(fields_structure{i})=Temp;
    percentages_test = (i/n_tests)*100;
    disp([num2str(percentages_test,2),'per cent completed'])
end
toc
bla=1.0;
end


function [Interpolation] = Interpolate_Data(Testdata,Testdata_a,Benchmark)
% 
% Input 
%==========================================================================
% Testdata = Dimensional results 
% Testdata_a = Dimensional results
%==========================================================================
% Output
% Interpolation vectors
%==========================================================================
    D_D = Testdata.D_norm;
    t_D = Testdata.time;
    D_A = Testdata_a.D_norm;
    t_A = Testdata_a.time; 
    D_Di = interp1(t_D,D_D,t_A);
    D_Ai = interp1(t_A,D_A,t_D);
    error_DA = nanmean(abs(D_D-D_Ai));
    error_AD = nanmean(abs(D_Di-D_A));
    if (Benchmark==1)
        disp([':::::::::::::::::::::::::::::::::::::::::::::::::::'])
        disp(['Errors interpolation are: ']);
        disp(['Adimensional interpolated to dimensional:',num2str(error_DA,'%3e')])
        disp(['Dimensional interpolated to adimensional:',num2str(error_AD, '%3e')])
        disp(['Error between detected detachment times:',num2str(abs(Testdata_a.t_det-Testdata.t_det),'%3e')])
        disp(['=================================================='])
    end
    Interpolation.t_A = t_A; 
    Interpolation.t_D = t_D;
    Interpolation.D_Di = D_Di; 
    Interpolation.D_Ai = D_Ai; 
    Interpolation.error_DA = [error_DA,nanmax(abs(D_D-D_Ai)),nanmin(abs(D_D-D_Ai))]; 
    Interpolation.error_AD = [error_AD,nanmax(abs(D_Di-D_A)),nanmin(abs(D_Di-D_A))]; 
end
function [Temp]=Processing_simulation(eta0,Df,n,l0,s0,D0,etaum,Benchmark,Df_UM,l_simulation)

    [ID] = Compute_slab_characteristics(eta0,Df,n,l0,s0,D0,etaum,Df_UM);
    % Run a simulation with a specific combination of parameter
    Testdata = Run_Simulation_Drag(ID,Benchmark);
    Testdata_a = Run_Simulation_DragA(ID.ID_A);
    %interpolation routine from AD->D
    [Interpolation] = Interpolate_Data(Testdata,Testdata_a,Benchmark);
    Temp = Testdata;
    Temp.Testdata_a=Testdata_a;
    Temp.initial_data = ID;
    Temp.Interpolation = Interpolation; 
    ID = [];
    Testdata = [];
end

