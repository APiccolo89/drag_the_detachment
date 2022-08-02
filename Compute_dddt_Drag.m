function Compute_dddt_Drag(D,B_n,B_d,n,L0,D0,drho,etaUM,tau0)
%Compute the D,t of the slab detachment. 
% Input Parameter: D0 = initial length 
%                  L0 = initial length
%                  drho = density contrast
%                  viscosity of the mantle
%                  B_d/n diffusion and dislocation creep compliance
%                  computed with the reference viscosity and stress 
% Alg: 
% => Compute the charcacteristic time
% => Compute the analytical solution detachment time 
% => Optimize the initial dDdt using fzero 
%                 a) Create a function handle to compute the initial
%                 guess(compute_guess_ddt)
%                 b) Create dDdt0 assuming no drag force
%                 c) Use Fzero and find dDdtA the initial dDdta 
% => Create a function handle to introduce into ode15i (variable method to
% solve non linear implicit function)
% => use and do plot
% =>save the associated vector in a structure and compute the next solution
%%
%Main Function

% Compute the refence strain rate 
ec = (B_n*tau0^n+B_d*tau0);
tc = 1/ec;
td = tc/n;
% Compute the initial guess of the dDdt 
% 1) Assuming that drag force is not active
dDdt0 = compute_guess_ddt(D,B_n,B_d,n,drho,L0,D0);
% 2) Assuming that drag force is active using optimizer
% compute the drag force
FUNFzero = @(x) compute_drag(x,D,B_n,B_d,n,etaUM,drho,D0,L0); 
% Compute the initial dDdt with drag force 
dDdtA = fzero(FUNFzero,dDdt0);
% Create function handle 
Funf_wi = @(t,x,xp0) compute_dragODE(x,xp0,B_n,B_d,n,etaUM,drho,D0,L0)
% Set the option for resolving the system of equation
options = odeset('RelTol',1e-2,'NormControl','on');
% resolve the system
[t,D] = ode15i(Funf_wi,[0 20*td],D0, dDdtA,options);

D_norm = D/D0;
tau_norm = 1./D_norm;

plot_time_evolution(t/td,D_norm,'D',append('D','_',num2str(L0),'_',num2str(D0),'_',num2str(etaUM)),n)
plot_time_evolution(t/td,tau_norm,'tau',append('tau','_',num2str(L0),'_',num2str(D0),'_',num2str(etaUM)),n)



end




function [dDdt0] = compute_guess_ddt(D,B_n,B_d,n,drho,L0,D0)

tau = (drho*D0*L0*9.81)/2/D;

dDdt0 = (-D*(B_n*tau^n+B_d*tau))/5;


end


function [res] = compute_drag(dDdt,D,B_n,B_d,n,etaUM,drho,D0,L0)

% For now the 1/h factor is assumed to be equal to the lenght of the slab


[tau_eff] = compute_the_effective_stress(D,dDdt,drho,D0,L0,etaUM,5.0);

res = -D*(B_n*(tau_eff)^n+B_d*(tau_eff))-dDdt; 


end

function [res] = compute_dragODE(D,dDdt,B_n,B_d,n,etaUM,drho,D0,L0)

% For now the 1/h factor is assumed to be equal to the lenght of the slab

[tau_eff] = compute_the_effective_stress(D,dDdt,drho,D0,L0,etaUM,5.0);

res = -D*(B_n*(tau_eff)^n+B_d*(tau_eff))-dDdt;
end


function plot_time_evolution(t,y,name,name_figure,n)
figure(1)
plot(t,y, "k",LineWidth=1.5);
grid on 
if strcmp(name,"D")
    hold on
    tt = [0:0.00001:1.0/n];
    anal = (1-n.*tt).^(1/n);
    plot(tt.*n,anal,'r')
    scatter(tt(400:400:end).*n,anal(400:400:end),'d')
end
xlim([0,4])
xlabel('t/td [n.d.]')
ylabel(append(name, '[n.d.]'))
print(name_figure,'-dpng')
clf; 
close; 
end

function [tau_eff] = compute_the_effective_stress(D,dDdt,drho,D0,L0,etaUM,alpha)
tau_B = (drho*D0*L0*9.81)/2/D;
tau_D = (4*etaUM*alpha*(D0^2/D^2)*dDdt)/2/D;
tau_eff = tau_B+tau_D;
end

