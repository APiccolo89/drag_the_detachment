clear all
close all

D0 = 80e3; 
L0 = 300e3; 
x  = 0:0.1e3:L0; 
gr = D0/L0; 
rhoS = 3367; 
rhoM = 3300; 
grho = rhoS/rhoM;
vd = -1e-14*D0; 
eta = 1e21; 
dP  = grho*gr*vd*(L0^2/(2*eta))^(-1) 
u = (1/(2*eta))*dP*(x.^2-x*L0)+vd/L0; 