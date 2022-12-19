
function [ID] = Compute_slab_characteristics(eta0DS,Df_S,n,l0,s0,D0,eta0DM,Df_UM,nlm)
% Output:
%==================================================================
% drho = density contrast
% B_d  = diffusion creep pre-exponential factor
% B_n  = dislocation creep pre-exponential factor
% tc   = characteristic time (s)
% ec   = characterstistc strain rate (1/s)
% Psi  = ratio between mantle effetive viscosity and slab effective
% viscosity
% Lambda = Combination of (l0/s)*(alpha*Psi)/2
%==================================================================
% Input
% eta0DS/M = reference viscosity at referenc stress
% s0       = reference stress
% l0       = initial lenght
% D0       = initial thickness
% n        = power law exponent
% Df_S/M   = viscosity contrast between diffusion and dislocation at reference
%           stress
% Description: from the primary input, the function retrieve the relevant
% field.It computes the dimensionless values of each parameter, and create
% additionally the structure through which the calculation are performed. 
%==================================================================
if nargin == 0 || nargin == 1 %default value for the unit test:
    eta0DS = 1.0e23;
    Df_S   = 10.0;
    n      = 3.5;
    l0     = 300e3;
    s0     = 100e6;
    D0     = 80e3;
    eta0DM  = 1e21;
    if nargin == 0
        nlm = Problem_type.NonLinear;
        Df_UM = 10.0;
    else
        nlm = Problem_type.Linear;
        Df_UM = 0.0;
    end
    % print the default value used for testing and all the other
    % parameters:
    disp ('Slab input parameter: ==================================')
    disp(['    D0 ( initial thickness thickness)   =' num2str(D0,1),' m' ])
    disp(['    L0 ( initial Length thickness)      =' num2str(l0,1),' m' ])
    disp(['    n ( stress exponent)                =' num2str(n,1),' [n.d.]' ])
    disp(['    s0 ( Bouyancy and reference stress) =' num2str(s0,1),' Pa' ])
    disp(['    eta0DS ( reference viscosity linear creep)    =' num2str(eta0DS,'%.E'),' Pas' ])
    disp(['    Df (viscosity contrast @ s0 )       =' num2str(Df_S,'%.E'),' [n.d.]' ])
    disp ('==============================================================')
    disp('UM input parameter:= =========================================')
    disp(['    etaum (  diffusion creep viscosity) =' num2str(eta0DM,'%.E'),' Pas' ])
    if nargin == 0
        disp(['    Df_UM (  viscosity contrast @s0 between dif/dis) =' num2str(Df_UM,'%.E'),' n.d.' ])
    else
        disp(['    Df_UM (  viscosity contrast @s0 between dif/dis) = is not set' ])
    end
    disp ('==============================================================')

end



eta0NS = (1/Df_S).*eta0DS; 
if nlm.islinear
    eta0NM = eta0DM; 
    n_m    = 1.0; 
else 
    eta0NM = (1/Df_UM).*eta0DM; 
end 
drho       = 2*s0/(9.81*l0);     % delta rho
B_n        = s0^(1-n)/eta0NS;       % compliance dislocation
B_d        = 1/(2.*eta0DS);         % compliance diffusion
ec         = (B_n*s0^n+B_d*s0);   % characteristic strain rate
tc         = 1/ec             ;   % characteristic time scale
etaS_eff   = (1/eta0NS+1/(eta0DS))^(-1); % effective viscosity of the slab at reference condition
Psi        = eta0DM/etaS_eff;    % ratio between the upper mantle viscosity and the slab viscosity
% Alpha
alpha      = 5.0;               % Ancient parameter derived by Yanick et al. 1986
s          = 1000e3;%1000e3*1.5;                    
len        = l0/(2*s);                     % Length divided by a characteristic lenght scale (i.e. size of my model)
Lambda     = len*alpha*Psi;             % Parameter derived by 2D numerical simulation
B_n_um     = (s0^(1-n))/(2*eta0NM);
ID.eta_CF = 1e18; % switch to 1.0 when the use wants to not deal with it. 
ID.B_D_C      = 1./(2.*ID.eta_CF);
ID.nlm  = nlm;

string_ID = {'B_d','B_n','s0','n','eta0DS','eta0NS','Df_S','drho','D0','l0','eta0DM','eta0NM','tc','ec','Psi','Lambda','alpha','len','Df_UM','s','B_n_um'};
if nargin == 0 || nargin == 1
    disp(['============================================================='])
    disp(['DIMENSIONAL       ==========================================='])
end
for is = 1:numel(string_ID)
    ID.(string_ID{is}) = eval(string_ID{is});
    str = strcat(string_ID{is}, ' = ', num2str(eval(string_ID{is})));
    if(nargin == 0 || nargin == 1)
        disp(['    ',str])
    end

end
ID.ID_A = Compute_slab_Adimensionals_IV(ID);

if nargin == 0 || nargin == 1
    disp(['ADIMENSIONAL      ==========================================='])

    fies = fields(ID.ID_A);
    for is = 1:numel(fies)
        str = strcat(fies{is}, ' = ', num2str(ID.ID_A.(fies{is})));
        disp(['    ',str])
    end
    disp(['============================================================='])
ID.ID_A.nlm       = nlm; 

end


end

function [ID_A] = Compute_slab_Adimensionals_IV(ID)
% Compute B_d,B_n,tauB,0,D_0, so forth adimensional
% B_nA = B_n/(tau{B,0}^-n*t_c^(-1))
% B_dA = B_d/(tau{B,0}^-1*t_c^(-1))
t_c = ID.tc; 
D0  = 1.0; 
tau_B0 = 1.0; 
B_d   = ID.B_d/((ID.s0*t_c)^(-1));
B_n   = ID.B_n/((ID.s0)^(-ID.n)*t_c^(-1)); 
B_dC  = ID.B_D_C/((ID.s0*t_c)^(-1));
Lambda = ID.Lambda;
% compute dDdt using dimensional value; 
dDdtBDC = ID.D0*(ID.B_d*ID.s0+ID.B_n*(ID.s0)^ID.n); 
dDdtADC = dDdtBDC*(1/(ID.D0))*(t_c);
%Compute the dDdt using adimensional value
dDdtA = 1.0*(B_d*1.0+B_n*(1.0)^ID.n);
if abs(dDdtA - dDdtADC) >1e-8 
    error('Issue with the non dimensionalization routine');
end
ID_A.alpha  = ID.alpha; 
ID_A.D0     = D0;
ID_A.s      = ID.s/ID.D0; 
ID_A.L0     = ID.l0/ID.D0; 
ID_A.B_d    = B_d;
ID_A.B_n    = B_n; 
ID_A.tau0   = tau_B0; 
ID_A.dDdtB  = dDdtA; 
ID_A.tc     = 1.0;
ID_A.n      = ID.n; 
ID_A.Lambda = Lambda; 
ID_A.etaS   = ((1/ID.eta0NS+1/(ID.eta0DS))^(-1))/(ID.s0*ID.tc);
ID_A.B_n_um = ID.B_n_um/((ID.s0)^(-ID.n)*t_c^(-1)); 
ID_A.eta0DM = ID.eta0DM/(ID.s0*ID.tc);
ID_A.eta_CF = ID.eta_CF/(ID.s0*ID.tc);
ID_A.Df_UM  = ID.Df_UM;
ID_A.Df_S   = ID.Df_S; 
ID_A.B_D_C = B_dC; 
ID_A.eta0DS = ID.eta0DS./(ID.s0*ID.tc);
ID_A.fetch(1)  = 0.0; 
ID_A.fetch(2)  = 1.0; 


end