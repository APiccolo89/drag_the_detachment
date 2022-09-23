function [ID] = Compute_slab_characteristics(eta0,Df,n,l0,s0,D0,etaum,Df_UM)
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
        % eta0 = reference viscosity at referenc stress
        % s0   = reference stress
        % l0   = initial lenght
        % D0   = initial thickness
        % n    = power law exponent
        % Df   = viscosity contrast between diffusion and dislocation at reference
        %        stress
        %==================================================================
        drho   = 2*s0/(9.81*l0);     % delta rho 
        B_n   = s0^(1-n)/eta0;       % compliance dislocation
        B_d   = 1/(Df*eta0);         % compliance diffusion
        ec    = (B_n*s0^n+B_d*s0);   % characteristic strain rate
        tc    = 1/ec             ;   % characteristic time scale
        etaS_eff = (1/eta0+1/(Df*eta0))^(-1); % effective viscosity of the slab at reference condition
        Psi      = etaum/etaS_eff;    % ratio between the upper mantle viscosity and the slab viscosity
        % Alpha 
        alpha = 5.0;                 % Ancient parameter derived by Yanick et al. 1986
        len = l0/(2*1000e3);         % Length divided by a characteristic lenght scale (i.e. size of my model)
        Lambda = len*alpha*Psi;      % Parameter derived by 2D numerical simulation 
        string_ID = {'B_d','B_n','s0','n','eta0','Df','drho','D0','l0','etaum','tc','ec','Psi','Lambda','alpha','len','Df_UM'};                    
        for is = 1:numel(string_ID)
              ID.(string_ID{is}) = eval(string_ID{is});
        end
        ID.eta_um_ = ID.etaum; 

        ID.ID_A = Compute_slab_Adimensionals_IV(ID);
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
Lambda = ID.Lambda;
% compute dDdt using dimensional value; 
dDdtBDC = ID.D0*(ID.B_d*ID.s0+ID.B_n*(ID.s0)^ID.n); 
dDdtADC = dDdtBDC*(1/(ID.D0))*(t_c);
%Compute the dDdt using adimensional value
dDdtA = 1.0*(B_d*1.0+B_n*(1.0)^ID.n);
if abs(dDdtA - dDdtADC) >1e-8 
    error('Issue with the non dimensionalization routine');
end
ID_A.D0 = D0;
ID_A.B_d = B_d;
ID_A.B_n = B_n; 
ID_A.tau0 = tau_B0; 
ID_A.dDdtB = dDdtA; 
ID_A.tc    = 1.0;
ID_A.n     = ID.n; 
ID_A.Lambda = Lambda; 
ID_A.Lambda_ = Lambda; % I need a buffer to remind it during the iteration
ID_A.Df_UM = ID.Df_UM;
end