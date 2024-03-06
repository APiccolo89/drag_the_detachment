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
    [Lambda,eta_um,eta_um_min,eta_um_MAX,tau_mantle,eps_mantle] = Compute_Lambda_re(ID_A,dDdt,D,nlm);    %Compute Lambda non linear
    tau_D = Lambda.*tau_B.^3.*dDdt.*ID_A.fetch(2);       % compute drag stress
end
tau_eff = tau_fetch.*(tau_B-(abs(tau_D)));


if nlm.islinear == 0
    varargout{1} = tau_eff;
    varargout{2} = tau_B;
    varargout{3} = tau_D;
    varargout{4} = Lambda;
    varargout{5} = eta_um;
    varargout{6} = eta_um_min;
    varargout{7} = eta_um_MAX;
    varargout{8} = tau_mantle;
    varargout{9} = eps_mantle;
    
    
else
    varargout{1} = tau_eff;
    varargout{2} = tau_B;
    varargout{3} = tau_D;
    varargout{4} = ID_A.Lambda;
end
end


function [Lambda,eta_um,eta_um_min,eta_um_MAX,tau,eps_mantle] = Compute_Lambda_re(ID_A,dDdt,D,nlm)
%==========================================================================
%Input: ID_A structure with the non dimensionalised parameter
%       dDdt non dimensionalised necking rate
%       D Actual thickness of the slab
%==========================================================================
%Output Scaled Lambda: Lambda_0*(1/(xium(gamma2*(tau_B^2)abs(dDdt))^(1/n-n)+1))
% => Divided the necessary variable, introduce an intermediate variable to
% avoid bug.
%==========================================================================

eta_um_min=[];
eta_um_MAX = [];
eps_mantle =[];
tau = [];



n      = ID_A.n;         %stress exponent
xium   = ID_A.Df_UM;     %xium
D0     = 1.0;            %D0
alpha  = ID_A.alpha;     %alpha
s      = ID_A.s;         %convective length scale
theta = (D0*alpha)./(s); %geometric parameter
n_cor    = (n-1)./n;
eps_tot=(theta.*(1./D).^2.*abs(dDdt));
temporary_var = 1+xium.*(eps_tot).^n_cor;
if nlm.iteration ==0
    eta_um    = ID_A.eta0DM./temporary_var;
    Lambda        = ID_A.Lambda./temporary_var;
    buffer_debug = (1+ID_A.Df_S).*(eta_um./ID_A.eta0DS).*(ID_A.L0*ID_A.alpha)./(2.*ID_A.s);
    
elseif nlm.iteration == -1
    
    eta_um_disl=1/2.*(ID_A.B_n_um).^(-1/n).*eps_tot.^(1/n-1);
    eta_um_diff  = ID_A.eta0DM;
    for i=1:length(eta_um_disl)
        eta_um_MAX(i)=min(eta_um_disl(i),eta_um_diff);
    end
    eta_um_min = ID_A.eta0DM./temporary_var;
    eta_um     = (eta_um_MAX+eta_um_min)./2;
    Lambda = (1+ID_A.Df_S).*(eta_um./ID_A.eta0DS).*(ID_A.L0*ID_A.alpha)./(2.*ID_A.s);
end


if nlm.iteration > 0
    for i = 1:length(eps_tot)
        eta_um_diff  = ID_A.eta0DM;
        eta_um_disl  = 1/2*(ID_A.B_n_um)^(-1/n)*eps_tot(i)^(1/n-1);
        eta_um_cut_off = 1/2*(1/2/ID_A.eta_CFU)^(-1)*eps_tot(i);
        eta_um_min(i)=(1/eta_um_diff+1/eta_um_disl)^(-1);
        eta_um_MAX(i) = min(eta_um_diff,eta_um_disl);
        tau_min = 2*eta_um_min(i)*eps_tot(i);
        tau_max = 2*eta_um_MAX(i)*eps_tot(i);
        
        tau_trial = [tau_min,tau_max];
        if abs(tau_trial(2)-tau_max(1))<1e-3
            tau_trial = (tau_trial(2)+tau_trial(1))/2;
        end
        options = optimset('TolX',1e-10);
        
        tau_f_h = @(x) compute_tau_mantle(ID_A,eps_tot(i),x,nlm);
        
        tau(i)     = fzero(tau_f_h,tau_trial,options);
     
        
        if nlm.cut_off ==1
            tau(i) = tau(i)+2*ID_A.eta_CF*eps_tot(i);
        end
        
        
        eta_um(i) = 0.5.*(tau(i)/eps_tot(i));
        if nlm.iteration == -1
            eta_um(i) =eta_um_MAX(i);
            eta_um_min=[];
            eta_um_MAX = [];
            
        else
            eta_um(i) = 0.5.*(tau(i)/eps_tot(i));
        end
        Lambda(i) = (1+ID_A.Df_S).*(eta_um(i)./ID_A.eta0DS)*(ID_A.L0*ID_A.alpha)/(2.*ID_A.s);
        eps_mantle(1,i)=(ID_A.B_n_um).*tau(i)^(n);
        eps_mantle(2,i)=(1/2/ID_A.eta0DM).*tau(i);
        eps_mantle(3,i)= eps_tot(i);
        eps_mantle(4,i)= eps_tot(i)-eps_mantle(1,i)-eps_mantle(2,i);
        %Lambda(i) = ID_A.Lambda./ID_A.eta0DM;
        %Lambda(i) = Lambda(i).*eta_um(i);
        
    end
    eta_um = eta_um';
    Lambda = Lambda';
else
    
end
end

% function
function [res] = compute_tau_mantle(ID,eps,tau,nlm)

B_n_m = ID.B_n_um;
B_d_m = 1./(2*ID.eta0DM);
n     = ID.n;

if nlm.cut_off == 0
    res = (eps-(B_d_m*(tau)+B_n_m*(tau)^n));
    
else
    B_d_u = 1/(2*ID.eta_CFU);
    res = (eps-(B_d_m*tau+B_n_m*tau^n+B_d_u*tau));
end
    res = res./eps;
end

