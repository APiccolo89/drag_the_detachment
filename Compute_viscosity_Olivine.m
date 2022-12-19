%=========================================================================%
% Depth dependent viscosity calculator. 
% ========================================================================%
% [INPUT]: -> [t0, age, Tp, D0,L0,(Vn,Vd)]=> Most difficult parameter     |
% 1st Construct a function that compute the reference viscosities for a   |
% given temperature of reference, and pressure of reference.              |
% 2nd Plot the reference viscosities corrected with the pressure          |
% [OUTPUT]:-> [eta0DM,eta0M,eta0DS,eta0S] = function                      |
%=========================================================================%
clear all
close all

% Reference values: 
Tp    = 1250+273.15; 
Pr    = 3300*100e3*9.81; 
t0    = 100e6;
Vnv   = [10e-6:1e-6:20e-6];
Vdv   = [5e-6];
age   = 10; 
D0 = 80e3;
L0 = [1,300,500,600].*1e3;
% Dry Olivine Data: 
[eta0DMP,eta0MP,eta0DSP,eta0SP] =  main_function(t0, age, Tp, Pr,D0,L0,Vnv,Vdv);


function [eta0DMP,eta0MP,eta0DSP,eta0SP] = main_function(t0, age, Tp,Pr, D0,L0,Vnv,Vdv)

[MMV,MML0]  = meshgrid(Vnv,L0);
eta0DMP = (MMV./MMV).*0.0;
eta0MP  = (MMV./MMV).*0.0;
eta0DSP = (MMV./MMV).*0.0;
eta0SP  = (MMV./MMV).*0.0;
xiumP   = (MMV./MMV).*0.0;
xium0   = (MMV./MMV).*0.0;

for i=1:length(MMV(:))
    l0 = MML0(i);
    Vn = MMV(i);

    [eta0DMP(i), eta0MP(i), eta0DSP(i), eta0SP(i),xiumP(i),xium0(i)] = compute_effective_viscosities(t0,D0,l0,Vn,Vdv,Tp,Pr,age);



end
bla = 0
figure(1)
contourf(MMV,MML0./1e3,log10(eta0DMP./eta0MP),5);colorbar; shading interp; colormap(crameri('nuuk',20));


figure(2)
La_0 = (MML0.*5.0)./(2*1e3);
xiS  = eta0DSP./eta0SP; 
Lambda_0 = ((1+xiS).*eta0DMP)./(eta0DSP);
Lambda_0(log10(Lambda_0)>0)=NaN; 
contourf(MMV,MML0./1e3,log10(Lambda_0),20);colorbar; shading interp; colormap(crameri('nuuk',20));
figure(3)
eta_eff_M  = 1./(1./eta0MP);%+1./eta0MP);
contourf(MMV,MML0./1e3,log10(eta_eff_M),20);colorbar; shading interp; colormap(crameri('nuuk',20));


end 
function [eta0MDP,eta0MP,eta0DSP,eta0SP,xium_P,xium_0] = compute_effective_viscosities(t0,D0,L0,Vn,Vd,Tp,Pr,age)  
 
% Constants
rho  = 3300; 
drho = (t0*2)/(9.81*L0);
alpha = 3e-5; 
Cp    = 1050; 
g = 9.81; 
R = 8.314;
phi = alpha/(Cp*rho);
%Dry Olivine Data
Bn = 6.22254e-16;
Bd = 1.5e-9; 
En = 530e3; 
n  = 3.5; 
Ed = 375e3;
Cn = (En*phi-Vn*(1-phi*Pr))/(R*Tp);
Cd = (Ed*phi-Vd*(1-phi*Pr))/(R*Tp);
% Slab length and thickness
 
% Compute Reference viscosity 
eta0M = compute_reference_viscosity(Bn,En,Vn,n,Tp,Pr,t0,R);
eta0DM = compute_reference_viscosity(Bd,Ed,Vd,1,Tp,Pr,t0,R);
xium_0 = eta0DM/eta0M; 
disp('=====================================================================')
disp(['Reference viscosity upper mantle at T = ' num2str(Tp), 'K, and at P = ',num2str(Pr/1e9,2), ' GPa and at tau0 = ' num2str(t0), 'MPa is:'])
disp(['log10(eta0) = ',num2str(log10(eta0M),4)])
disp(['log10(eta0D) = ', num2str(log10(eta0DM),4)])
disp('=====================================================================')
% Compute Reference viscosity slab 
Tavg =  compute_average_T_slab(D0,age,Tp);
eta0S = compute_reference_viscosity(Bn,En,Vn,n,Tavg,Pr,t0,R);
eta0DS = compute_reference_viscosity(Bd,Ed,Vd,1,Tavg,Pr,t0,R);
xiuS_0 = eta0DS/eta0S; 
% disp('=====================================================================')
% disp(['Reference viscosity Slab at T = ' num2str(Tavg), 'K, and at P = ',num2str(Pr/1e9,2), ' GPa and at tau0 = ' num2str(t0), 'MPa is:'])
% disp(['log10(eta0) = ',num2str(log10(eta0S),4)])
% disp(['log10(eta0D) = ', num2str(log10(eta0DS),4)])
% disp('=====================================================================')
% Compute integral Upper mantle 
Exn = @(x) compute_exponential(Cn,x,phi,3300*g);
Exd = @(x) compute_exponential(Cd,x,phi,3300*g);
intexn = (integral(Exn,0,L0))./L0;
intexd = (integral(Exd,0,L0))./L0;
eta0MP = eta0M*intexn;
eta0MDP = eta0DM*intexd; 
xium_P = eta0MDP/eta0MP; 
if isnan(xium_P)
bla = 0; 
end
disp('=====================================================================')
disp(['Average reference viscosity UM integrated along the slab = '])
disp(['log10(eta0) = ',num2str(log10(eta0MP),4)])
disp(['log10(eta0D) = ', num2str(log10(eta0MDP),4)])
disp(['xiUM0 = ' num2str(log10(xium_0)), ' and xiUMP =', num2str(log10(xium_P))])
disp('=====================================================================')
% Compute integral Slab 
CnS = (En*phi-Vn*(1-phi*Pr))/(R*Tavg);
CdS = (Ed*phi-Vd*(1-phi*Pr))/(R*Tavg);
ExnS = @(x) compute_exponential(CnS,x,phi,(rho+drho)*g);
ExdS = @(x) compute_exponential(CdS,x,phi,(rho+drho)*g);
intexnS = (integral(ExnS,0,L0))./L0;
intexdS = (integral(ExdS,0,L0))./L0;
eta0SP = eta0S*intexn;
eta0DSP = eta0DS*intexd; 
xiS_P = eta0DSP/eta0SP; 
% Compute the viscosity contrast 
end
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

function [exp_] = compute_exponential(C,l,phi,w)
    exp_= exp(-C.*((l.*w)./(1+phi.*l*w)));
end

