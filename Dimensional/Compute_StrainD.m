function [eps_eff,eps_dif,eps_dis] = Compute_StrainD(ID,tau_eff,Benchmark)
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
    tau_eff = 100e6;
    Benchmark = 0; 
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

[eps_eff] = compute_effective_viscosity(tau_eff,ID,eps_eff);


%Benchmark 
% Small function that call the adimensional computation of strain rate and
% uses the data coming from the computation of tau_eff and compare with
% actual data

if  Benchmark == 1

    [A,B,C]=Compute_StrainA(ID.ID_A,tau_eff./ID.s0); 
    res_DimB = sum(eps_eff./ID.ec-A);
    res_Dimn = sum(eps_dis/ID.ec-C);
    res_Dimd = sum(eps_dif/ID.ec-B);
    if length(A)>1
        disp([':::::::::::::::::::::::::::::::::::::::::::::::::::'])
        disp(['Dimensional-Adimensional computation of effective strain has an error of:'])
        disp(['effective strain is ',num2str(res_DimB,'%10.5e')])
        disp(['dislocation creep  is ',num2str(res_Dimn,'%10.5e')])
        disp(['diffusion creep  is ',num2str(res_Dimd,'%10.5e')])
        disp([':::::::::::::::::::::::::::::::::::::::::::::::::::'])

    end
    if (abs(res_DimB)>1e-6 )
      error('Error between the computations is unsustainable')
    end


end




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
function [eps_eff] = compute_effective_viscosity(tau_eff,ID,eps_eff)

eta_eff = ID.eta0DS./(1+ID.Df_S.*tau_eff.^(ID.n-1));
if eta_eff < ID.eta_CF
    eps_eff = ID.B_D_C.*tau_eff; 
else 
    eps_eff = eps_eff; 
end

end

