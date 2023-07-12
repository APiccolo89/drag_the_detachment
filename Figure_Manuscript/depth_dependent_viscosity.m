close all;
clf; 
clear all; 

Tp = 1300+273.15; 
Tp4 = 1000+273.15; 
Tp3 = 800+273.15;
Tp2 = 600+273.15; 
Pr   = 3200*9.81*100e3; 
tau0 = 100e6; 
L0   = 600e3; 
rho  = 3300; 
alpha = 3e-5; 
Cp    = 1050; 
g = 9.81; 
R = 8.314;
phi = alpha/(Cp*rho);
%==========================================================================
% Rheology 
%==========================================================================
Bn = 6.22254e-16;
Bd = 1.5e-9; 
En = 530e3; 
Vn = 15e-6; 
n  = 3.5; 
Ed = 375e3;
Vd = 5e-6; 
% Viscosity @ Reference P,T,tau
eta0 = (tau0^(1-n)/(2*Bn))*exp((En+Pr*Vn)/(R*Tp)); 
eta0D = (tau0^(1-1)/(2*Bd))*exp((Ed+Pr*Vd)/(R*Tp)); 
xium  = eta0D/eta0;
eta02 = (tau0^(1-n)/(2*Bn))*exp((En+Pr*Vn)/(R*Tp2)); 
eta0D2 = (tau0^(1-1)/(2*Bd))*exp((Ed+Pr*Vd)/(R*Tp2)); 
xiS2  = eta0D2/eta02;
eta03 = (tau0^(1-n)/(2*Bn))*exp((En+Pr*Vn)/(R*Tp3)); 
eta0D3 = (tau0^(1-1)/(2*Bd))*exp((Ed+Pr*Vd)/(R*Tp3)); 
xiS3  = eta0D3/eta03;
eta04 = (tau0^(1-n)/(2*Bn))*exp((En+Pr*Vn)/(R*Tp4)); 
eta0D4 = (tau0^(1-1)/(2*Bd))*exp((Ed+Pr*Vd)/(R*Tp4)); 
xiS4  = eta0D4/eta04;


Ll = 200e3:10e3:600e3; 
tdum = 0.0001:0.0001:1.0; 
gamma = (Ll.*5.0)./2000e3;
Psi2 = (1+xiS2)*eta0D/eta0D2;
Psi3 = (1+xiS3)*eta0D/eta0D3;
Psi4 = (1+xiS4)*eta0D/eta0D4;
for i=1:length(tdum)
    for j=1:length(gamma)
    Lambda2(i,j) =  (gamma(j).*Psi2)./(1+xium*(tdum(i))^(n-1));
    Lambda3(i,j) =  (gamma(j).*Psi3)./(1+xium*(tdum(i))^(n-1));
    Lambda4(i,j) =  (gamma(j).*Psi4)./(1+xium*(tdum(i))^(n-1));
    end
end

%Lambda2 = (gamma.*Psi2)./(1+xium*(0.001)^(n-1)); 
%Lamdba3 = (gamma.*Psi3)./(1+xium*(0.001)^(n-1)); 
%Lambda4 = (gamma.*Psi4)./(1+xium*(0.001)^(n-1));

figure(1)
contourf(gamma,log10(tdum.^(n-1)),log10(Lambda4),10);shading interp
colormap("jet")
caxis([-14 0])
colorbar; 
xlabel('\gamma')
ylabel('log_{10}(\tau^{n-1})')
figure(2)
plot(Ll./1e3,log10(Lambda2(1,:)))
hold on
plot(Ll./1e3,log10(Lambda3(1,:)))
plot(Ll./1e3,log10(Lambda4(1,:)))
%==========================================================================
% Main Vector
%==========================================================================
tau  = (0.01:0.1:100).*1e6; 
TT    = tau./tau0; 
l    = 0:100:L0;  
T = Tp*(1+phi*l*g*rho); 
P = Pr+l*g*rho;
Cn = (En*phi-Vn*(1-phi*Pr))/(R*Tp);
Cd = (Ed*phi-Vd*(1-phi*Pr))/(R*Tp);
detan = zeros(length(P),length(TT));
%==========================================================================
% Function viscosities
%==========================================================================
% Compute the viscosity profile diffusion and dislocation as a function of
% the rheological properties and givens deviatoric stress tensor second
% invariant. 
%-------------------------------------------------------------------------%
exp_D = compute_exponential(Cd,l,phi,rho*g); 
exp_n = compute_exponential(Cn,l,phi,rho*g); 
etaD  = visc_calc(eta0D,1,exp_D,tau0);  
eta_n(1,:) = visc_calc(eta0,n,exp_n,1e6/tau0);
eta_n(2,:) = visc_calc(eta0,n,exp_n,10e6/tau0);
eta_n(3,:) = visc_calc(eta0,n,exp_n,50e6/tau0);
eta_n(4,:) = visc_calc(eta0,n,exp_n,100e6/tau0);
eta_eff(1,:)= (1./eta_n(1,:)+1./etaD).^(-1);
eta_eff(2,:)= (1./eta_n(2,:)+1./etaD).^(-1);
eta_eff(3,:)= (1./eta_n(3,:)+1./etaD).^(-1);
eta_eff(4,:)= (1./eta_n(4,:)+1./etaD).^(-1);
%-------------------------------------------------------------------------%
% Plot simplified picture 
%-------------------------------------------------------------------------%
figure(1)
plot((etaD),-l./1e3,LineWidth=1.5,Color='red',LineStyle='--')
hold on
plot(eta_n(1,:),-l./1e3,LineWidth=0.8,Color='k',LineStyle=':')
plot(eta_n(2,:),-l./1e3,LineWidth=0.8,Color='b',LineStyle=':')
plot(eta_n(3,:),-l./1e3,LineWidth=0.8,Color='g',LineStyle=':')
plot(eta_n(4,:),-l./1e3,LineWidth=0.8,Color='c',LineStyle=':')





plot(eta_eff(1,:),-l./1e3,LineWidth=1.5,Color='k',LineStyle='-')
plot(eta_eff(2,:),-l./1e3,LineWidth=1.5,Color='b',LineStyle='-')
plot(eta_eff(3,:),-l./1e3,LineWidth=1.5,Color='g',LineStyle='-')
plot(eta_eff(4,:),-l./1e3,LineWidth=1.5,Color='c',LineStyle='-')

legend('\eta_D, [Pas]','\eta_n \tau = 1 MPa, [Pas]','\eta_n \tau = 10 MPa, [Pas]','\eta_n \tau = 50 MPa, [Pas]','\eta_n \tau = 100 MPa, [Pas]',lcn='southwestoutside')

xlabel("$log_{10}(\eta_{eff}) [Pas]$",Interpreter="latex")
ylabel("$l_{slab} [km]$",Interpreter="latex")
set(gca, 'XScale', 'log')
grid on 
%-------------------------------------------------------------------------%
%=========================================================================%
%Compute integrals of exponential
%=========================================================================%


int(1) =  integral_exp(0,rho*g,Cd,phi); 
int(2) =  integral_exp(l(end),rho*g,Cd,phi);
V = int(2)-int(1); 
eta_M  = (eta0D*V)/300e3; 


detaD = eta0D.*exp(-Cd.*((l.*g*rho)./(1+phi.*l*g*rho)));

[M] = compute_numerical_integral(detaD,l);

for i = 1:length(TT)
    detan(:,i) = (eta0).*(TT(i))^(1-n).*exp(-Cn.*((l.*g*rho)./(1+phi.*l*g*rho)));
    ratio(:,i) = detan(:,i)./detaD'; 
    eff(:,i)   = (1./detan(:,i)+1./detaD').^-1;
end
etaS  =  (tau0^(1-n)/(2*Bn))*exp((En+Pr*Vn)/(R*1073));
dDdt = 8e-10; 
tau_D= (dDdt/300e3)*5.0.*eff;



%eta_l = eta0.*exp(-C.*((l.*g*rho)./(1+phi.*l*g*rho))); 
%eta_P = eta0.*exp((E./(R.*T).*(1-(T./Tp)))).*exp(V./(R.*T).*(P-Pr.*(T./Tp)));




figure(1)
plot(l/1e3,log10(detan(:,1:1:end)),'Color','red')
hold on
plot(l/1e3,log10(detaD),'Color','blue',LineStyle=':',LineWidth=3.0)

figure(3)
plot(l/1e3,log10(eff(:,1:1:end)),'Color','red')
hold on
plot(l/1e3,log10(detaD),'Color','blue',LineStyle=':',LineWidth=3.0)

level=0.1:1:6;
figure(2)
contourf(l./1e3,tau./tau0,log10(eff'),20)
hold on
%plot(l./1e3,tau_D'./tau0,LineStyle=":",Color='k',LineWidth=1.2)
ylim([0,1])
colormap("jet")
colorbar


figure(3)
contourf(tau./tau0,log10(eff'),tau_D,20)
hold on
%plot(l./1e3,tau_D'./tau0,LineStyle=":",Color='k',LineWidth=1.2)
colormap("jet")
colorbar

function [exp_] = compute_exponential(C,l,phi,w)
    exp_= exp(-C.*((l.*w)./(1+phi.*l*w)));
end


function [eta] = visc_calc(eta0,n,exp,tau) 
    eta = eta0.*tau^(1-n).*exp; 
end


function [V] = integral_exp(l0,w,C,phi)
V = (exp(-C)*((l0*phi+1)*exp(C/(l0*phi + 1))-C*expint(C/(l0*phi+1))))/phi; 
end

function [M] = compute_numerical_integral(eta,l)
dl = diff(l);
eta05 = 0.5.*(eta(2:1:end)+eta(1:1:end-1))
M   = sum(eta05.*dl(1))/(l(end)-l(1)); 


end



