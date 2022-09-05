% compute the change in interface roughness
%----------------------------------------------------
function drdt = Compute_drdt(Sol,D,r_sdis0,deltaQg,deltaQdif,n,m,q,T0)

% input:
% Sol: vector, consisting of:
%  - interface curvature r
%  - temperature perturbation T
%  - slab thickness d
%
% nondimensional material parameters: 
% D, r_sdis0, deltaQg,deltaQdif, n,m,T0


r = Sol(1);
T = Sol(2);
d = Sol(3);

% precompute inverse of thickness, as this is mostly used
invd = 1./d;

% precompute temperature factor
Tfac = T./(1+T./T0);

% term1: growth
term1 =  r_sdis0.^(q+1) .* exp(deltaQg*Tfac) ./r.^(q-1);

% term2: reduction due to dislocation creep
term2 = exp(Tfac).* invd.^(n+1) .* r.^2;

% term3: reduction due to diffusion creep
term3 = exp(deltaQdif.*Tfac).*invd.^2 .* r.^(2-m);

% everything together
drdt = D.*(term1-term2-term3);
