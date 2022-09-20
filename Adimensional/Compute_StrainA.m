function [eps_eff,eps_dif,eps_dis] = Compute_StrainA(ID,tau_eff,Benchmark)
% Input: 
%==========================================================================
% ID, data structure with all the initial data, both
% dimensional/adimensional (should introduce a flag, because the
% computation of the strain rate should be the same in both the approach)
% tau_eff is the effective stress compute
%==========================================================================
% Output: 
% eps_eff effective strain rate (eps_dif+eps_dis) and the two main
%==========================================================================
% component of the effective strain rate. 
% Function => Unit test: computed B_d,B_n in excel, retrieve the
% fundamental data => strain rate should be 1.0 in adimensional. 
%==========================================================================

if nargin == 0 % Unit test
    %Default value compute using 100e6 reference stress, and reference
    %viscosity 10e21 dislocation, 10e22 diffusion creep. 
    B_n = 5.00E-42;
    B_d = 5.00E-23;
    n   = 3.5;
    check_eff = 5.50E-14;
    %non dimensional 
    B_d = B_d/((100e6)^(-1)*check_eff);
    B_n = B_n/((100e6)^(-n)*check_eff);
    check_eff=check_eff/check_eff;
    tau_eff = 1.0; 
else
    B_n = ID.B_n;
    B_d = ID.B_d;
    n   = ID.n  ;
end
% peirls creep place holder
%B_p = ID.B_p; 

% Compute dislocation creep 
eps_dis = Compute_strain_m(B_n,tau_eff,n); 
% Compute diffusion creep 
eps_dif = Compute_strain_m(B_d,tau_eff);
% Compute peirls creep 
%place holderai

if nargin == 0 
    if abs(eps_dif+eps_dis-check_eff)>1e-5
        error('Something wrong with Compute_StrainA')
    else
        disp('Nothing to declare, I am a clean function')
    end

end
% Compute the effective strain 
eps_eff = eps_dif+eps_dis; 
end

function [eps] = Compute_strain_m(B,tau_eff,n)
% For a given B, tau_eff, and n, gives the strain rate 
%==========================================================================
% Input: B Pre-exponential factor 
%        n stress exponent
%        tau_eff the effective stress 
%        if n is not one of the output parameter, default = 1.0 =>
%        Diffusion creep
%==========================================================================
% Output: epsilon_rate for a given mechanism
%==========================================================================
if nargin == 2
    n = 1.0;
end
eps = B.*tau_eff.^n;    
end
