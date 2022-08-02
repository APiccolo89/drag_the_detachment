%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clf



eta0 = 1e23; 
eta0D = 10*eta0;
L0    = 500e3;
D0    = 100e3; 
n     = 3.5; 
tau0  = 100e6;
drho   = 2*tau0/(9.81*L0);
B_n   = tau0^(1-n)/eta0;
B_d   = 1/eta0D; 
etaUM = 1e22; 
name_1 = num2str(log10(eta0));

D = D0;
Compute_dddt_Drag(D,B_n,B_d,n,L0,D0,drho,etaUM,tau0,name_1); 