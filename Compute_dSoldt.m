% compute the solution change for one time/iteration step

function dSoldt = Compute_dSoldt(t,Sol,D,S,r_sdis0,deltaQg,deltaQdif,n,m,q,T0,tscale)

% change in interface roughness
drdt = Compute_drdt(Sol,D,r_sdis0,deltaQg,deltaQdif,n,m,q,T0);
% change in temperature perturbation
dTdt = Compute_dTdt(Sol,S,T0,deltaQdif,n,m);
% change in slab thickness
dddt = Compute_dddt(Sol,T0,deltaQdif,n,m);

% solution change
dSoldt = [drdt;dTdt;dddt];

% as changes are sometimes too large to be handled accurately, time is additionally
% scaled
dSoldt = dSoldt/tscale;