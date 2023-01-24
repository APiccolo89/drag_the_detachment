function [eta0] = compute_reference_viscosity(B,E,V,n,Tp,Pr,t0,R)
% Compute the reference viscosity for a given mechanism. Diffusion creep n
% ==1, while power law != 1. 
% =========================================================================
% Input: 
%==========================================================================
% B,E,V,n => Rheological law parameter (Pre-exponential factor, energy of
% activaction, volume of activation, stress exponent) 
%==========================================================================
% Output: 
%==========================================================================
% eta0 : reference viscosity for a given mechanism (which is going to be
% multiplied for the integral of the exponential and the normalized stress
% to obtain the actual viscosity
%==========================================================================
eta0 = 0.5*(1/B)*t0^(1-n)*exp((E+Pr*V)/(R*Tp));
end
function [avgT] = compute_average_T_slab(D0,age,Tp)
kappa = 1e-6; 
D0v    = 0:0.5e3:D0; 
age = age*(365.25*60*60*24*1e6); 
T    =  half_space(D0v,age,Tp,kappa);
fun  = @(x) half_space(x,age,Tp,kappa);
int = (integral(fun,0,80e3))./80e3; 
avgT = int; 
end
function [T] = half_space(D0v,age,Tp,kappa)
T = Tp-Tp.*erfc(D0v./2./sqrt(kappa*age));
end


