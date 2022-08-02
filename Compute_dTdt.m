% compute the change in temperature perturbation
%----------------------------------------------------
function dTdt = Compute_dTdt(Sol,S,T0,deltaQdif,n,m)

% input:
% Sol: vector, consisting of:
%  - interface curvature r
%  - temperature perturbation T
%  - slab thickness d
%
% nondimensional material parameters: 
% S,,deltaQdif, n,m,T0

r = Sol(1);
T = Sol(2);
d = Sol(3);

% precompute inverse of thickness, as this is mostly used
invd = 1./d;

% precompute temperature factor
Tfac = T./(1+T./T0);

% term1: shear heating due to dislocation creep
term1 = exp(Tfac).* invd.^(n+1);

% term2: shear heating due to diffusion creep
term2 = exp(deltaQdif.*Tfac).*invd.^2 .* r.^(-m);

% everything together
dTdt = S.* (term1 + term2);

