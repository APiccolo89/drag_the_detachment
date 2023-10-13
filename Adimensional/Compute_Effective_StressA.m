function varargout = Compute_Effective_StressA(D,dDdt,ID_A,nlm)
    % Compute the adimensional effective stress. 
    % Input ID_A: initial input adimensionalized with the characteristic
    % length scale, and computation of Lambda parameter that suits the
    % simulation with realistic parameter
    %=====================================================================
    % Eq. tau_eff = (F_eff)/2D => (F_B+F_D)/2D =>
    % F_B(1+F_D/F_B)=>tau_{0,b}*D0/D(1+Lambda(D0/D)^2*dDdt/dDdt|0,b)=> a.d. =>
    % tau_eff/tau_{0,b} => D0/D(1+Lambda(D0/D)^2*dDdt') D0/D = tau_B (a.d.)
    % => tau_B = D0/D (tau_B=F_B/2D => F_B=2tau_{B,0}D_0=> tau_B(a.d.) =
    % tau_B/tau_{B_0} => D0/D
    %======================================================================
    % out parameter => tau_eff, tau_B,tau_D adimensional
    % Buoyancy stress computed tau_B = F_B/2/D;
    % unit_test @ assuming Lambda = 0.001 D=D_0 and dDdt =0.9 what is the
    % stress
    %======================================================================

    tau_B = (ID_A.D0./D); % Buoyancy force
    tau_fetch = (1-ID_A.fetch(1)); % Depth correction
    % Drag force related stress compute formulation of Bercovici et al 2015 

    if nlm.islinear
        tau_D = +ID_A.Lambda.*tau_B.^3.*ID_A.fetch(2).*dDdt;  % Compute drag stress linear
    else
        %compute alternative lambda using dDdt 
        [Lambda,eta_um] = Compute_Lambda_re(ID_A,dDdt,D,nlm);    %Compute Lambda non linear
        tau_D = Lambda.*tau_B.^3.*dDdt.*ID_A.fetch(2);       % compute drag stress
    end
    tau_eff = tau_fetch.*(tau_B-(abs(tau_D)));

    if nlm.islinear == 0 
        varargout{1} = tau_eff;
        varargout{2} = tau_B; 
        varargout{3} = tau_D; 
        varargout{4} = Lambda; 
        varargout{5} = eta_um; 
    else
        varargout{1} = tau_eff;
        varargout{2} = tau_B; 
        varargout{3} = tau_D; 
        varargout{4} = ID_A.Lambda; 
    end
end


function [Lambda,eta_um] = Compute_Lambda_re(ID_A,dDdt,D,nlm)
%==========================================================================
%Input: ID_A structure with the non dimensionalised parameter
%       dDdt non dimensionalised necking rate
%       D Actual thickness of the slab
%==========================================================================
%Output Scaled Lambda: Lambda_0*(1/(xium(gamma2*(tau_B^2)abs(dDdt))^(1/n-n)+1))
% => Divided the necessary variable, introduce an intermediate variable to
% avoid bug. 
%==========================================================================
n      = ID_A.n;         %stress exponent
xium   = ID_A.Df_UM;     %xium
 D0     = 1.0;            %D0
 alpha  = ID_A.alpha;     %alpha
 s      = ID_A.s;         %convective length scale
 theta = (D0*alpha)./(s); %geometric parameter
 n_cor    = (n-1)./n; 
 eps_tot=(theta.*(D0./D).^2.*abs(dDdt)); 
 temporary_var = 1+xium.*(theta.*(D0./D).^2.*abs(dDdt)).^n_cor;
 eta_um    = ID_A.eta0DM./temporary_var; 
 if nlm.iteration==1
    tau_trial = 2*eps_tot*eta_um;




 end







if eta_um < ID_A.eta_CF & ID_A.cut_off_Mantle >0.0
    Lambda     = ID_A.Lambda./ID_A.eta0DM;
    Lambda     = Lambda.*ID_A.eta_CF; 
    eta_um    = ID_A.eta_CF;
else
    Lambda        = ID_A.Lambda./temporary_var; 
end
end

% function 
function [res] = compute_tau_mantle(ID,eps,tau)

    B_n_m = ID.B_n_um;
    B_d_m = 1./(2*ID.eta0DM);
    n     = ID.n; 
    if nlm.cut_off == 0
        res = eps-(B_d_m*tau+B_n_m*tau^n);
    else
        B_d_u = ID_A.B_D_upper; 
        res = eps-(B_d_m*tau+B_n_m*tau^n+ID_A.B_D_upper);
    end
end

