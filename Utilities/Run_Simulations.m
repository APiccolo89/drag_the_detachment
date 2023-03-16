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
disp(['=============================================================='])

disp(['Number of tests is ',num2str(n_tests)])
disp(['=============================================================='])

pause(3)
for i=1:n_tests
    T_name   = strcat('T_',num2str(i));
    eta0DM   = ium(i);
    s0       = iut0(i);
    eta0DS   = iuet0(i);
    l0       = iuL0(i);
    xiUS     = iuDfS(i);
    xiUM     = iuDfM(i);
    D0       = iuD0(i);
  %  try
        [Temp]   = Processing_simulation(eta0DS,xiUS,n,l0,s0,D0,eta0DM,Benchmark,xiUM,nlm);
   
%catch 
 %       disp([num2str(i),'out of',num2str(ntests), 'tests'])
  
%end
    Tests.(T_name)=Temp;
    
    disp([num2str(i),'out of',num2str(n_tests), 'tests'])
end
toc
end


