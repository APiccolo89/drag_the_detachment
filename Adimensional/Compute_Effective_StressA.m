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
    if nargin == 0
        ID_A.Lambda     = 0.001;
        ID_A.D0         = 1.0; 
        D               = 0.5;
        dDdt            = -1.5; 
        check_eff       = 1.988;
        check_drag      = -0.012;
        check_B         = 2.0; 
        nlm = Problem_type.Linear;
    end
   

    tau_B = (ID_A.D0./D);
    tau_fetch = (1-ID_A.fetch(1));%.*tau_B; 
    % Drag force related stress compute formulation of Bercovici et al 2015 

    if nlm.islinear
        tau_D = +ID_A.Lambda.*tau_B.^3.*ID_A.fetch(2).*dDdt; 
    else
        [tau_M,Lambda] = compute_drag_stressA(ID_A,D,dDdt); 
        tau_D = Lambda.*tau_B.^3.*dDdt.*ID_A.fetch(2);
    end
    tau_eff = tau_B.*tau_fetch+tau_D;

    % Effective stress
    if nargin == 0
        if abs(tau_eff-check_eff)>1e-5 || abs(tau_D-check_drag)>1e-5 || abs(tau_B-check_B)>1e-5
            error('Something wrong with function Compute_Effective_StressA')
        else
            disp('Nothing to declear, I am innocent')
        end
    end
    if nlm.islinear == 0 
        varargout{1} = tau_eff;
        varargout{2} = tau_B; 
        varargout{3} = tau_D; 
        varargout{4} = Lambda; 
    else
        varargout{1} = tau_eff;
        varargout{2} = tau_B; 
        varargout{3} = tau_D; 
        varargout{4} = ID_A.Lambda; 
    end
end
