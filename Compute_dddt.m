% compute the change in slab thickness
%----------------------------------------------------
function dddt = Compute_dddt(Sol,T0,deltaQdif,n,m)

% input:
% Sol: vector, consisting of:
%  - interface curvature r
%  - temperature perturbation T
%  - slab thickness d
%
% nondimensional material parameters: 
% deltaQdif, n,m,T0

r = Sol(1);
T = Sol(2);
d = Sol(3);

% precompute inverse of thickness, as this is mostly used
invd = 1./d;

% precompute temperature factor
Tfac = T./(1+T./T0);

% term1: shear heating due to dislocation creep
term1 = exp(Tfac).* invd.^(n-1);

% term2: shear heating due to diffusion creep
term2 = exp(deltaQdif.*Tfac).* r.^(-m);

% everything together
dddt = -1.* (term1 + term2);