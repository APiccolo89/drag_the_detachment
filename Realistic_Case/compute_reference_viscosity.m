function [eta0,D0v,eta_z,T_z] = compute_reference_viscosity(B,E,V,n,Tp,Pr,t0,R,age)
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
D0v=[];
eta_z=[];
T_z=[]; 


eta0 = 0.5*(1/B)*t0^(1-n)*exp((E+Pr*V)/(R*Tp));
if ~isnan(age)
    kappa = 1e-6;
    age = age*(365.25*60*60*24*1e6);
    fun  = @(x) half_space(x,age,1350,kappa);
    eta_ = @(x) eta_fun(x,B,E,V,n,Pr,t0,R,fun);  
    
    D0v    = 0:0.5e3:80e3;
    T_z    =  half_space(D0v,age,1350+273.15,kappa);
    eta_z = eta_(D0v); 
    int = (integral(eta_,0,80e3))./80e3;
    eta0 = int; 


end



end

function [T] = half_space(D0v,age,Tp,kappa)
T = Tp-(Tp-(20+273.15)).*erfc(D0v./2./sqrt(kappa*age));
end



function [eta0] = eta_fun(d,B,E,V,n,Pr,t0,R,fun) 
T = fun(d); 
eta0 = 0.5.*(1./B).*t0^(1-n).*exp((E+Pr.*V)./(R.*T));
end