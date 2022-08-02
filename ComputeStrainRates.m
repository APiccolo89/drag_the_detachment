% compute individual strain rates and contributions from grain size and
% temperature
function [edis,edis_T,edif,edif_T,edif_r] = ComputeStrainRates(r,T,d,deltaQdif,n,m,T0)

% precompute temperature factor
Tfac = T./(1+T./T0);
% precompute inverse of thickness, as this is mostly used
invd = 1./d;

% term1: shear heating due to dislocation creep
edis_T = exp(Tfac);
edis   = edis_T.* invd.^(n);

% term2: shear heating due to diffusion creep
edif_T = exp(deltaQdif.*Tfac);
edif_r = r.^(-m);
edif   = edif_T.*invd .* edif_r;