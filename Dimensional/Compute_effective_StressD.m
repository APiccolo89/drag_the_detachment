function varargout = Compute_effective_StressD(D,dDdt,ID,Benchmark,non_linear_um)
    % Compute effective stress
    % Input data: 
    % D => Actual thickness 
    % dDdt => actual necking rate
    % ID   => actual initial data
    %======================================================================
    % Output data: 
    % tau_eff effective stress
    % tau_B   buoyancy stress
    % tau_D   drag force
    %======================================================================
    if nargin == 0
        % Checking value are computed using Excell
        ID.D0 = 80e3;
        ID.s0 = 100e6; 
        D       = 40e3; 
        dDdt    = -ID.D0*1e-14; 
        ID.etaum = 1e20; 
        ID.len   = 300e3/1000e3;
        ID.alpha = 5.0;
        check_D  = -1.20E+07;
        check_B  = 2.00E+08;
        check_Eff = 1.88E+08;
        Benchmark = 0; 
    end
    % Buoyancy stress computed tau_B = F_B/2/D; 
    tau_B = (ID.s0*ID.D0)./D;
    % Drag force related stress compute formulation of Bercovici et al 2015 
    % Compute the effective stress and effective viscosity 
    
    if isnan(ID.Df_UM) || ~ID.Df_UM 
        tau_D = (2*ID.etaum*ID.alpha*(ID.D0^2./D.^2).*dDdt.*ID.len)./2./D;
    else
        [tau_D,eta_um,ID] = compute_drag_stress(ID,D,dDdt); 
    end
    % benchmark: 
    
    if  Benchmark == 1 
        
        [t_effA,t_BA,t_DA]=Compute_Effective_StressA(D/ID.D0,dDdt/(ID.D0/ID.tc),ID.ID_A,0.0); 
        res_DimB = sum(tau_B./ID.s0-t_BA);
        res_DimD = sum(tau_D./ID.s0-t_DA);

        if length(t_BA)>1
            disp(['=================================================='])
            disp(['Dimensional-Adimensional computation of effective stress has an error of:'])
            disp(['Bouyancy stress is ',num2str(res_DimB,'%10.5e')])
            disp(['Drag force is ',num2str(res_DimD,'%10.5e')])
            disp([':::::::::::::::::::::::::::::::::::::::::::::::::::'])
        end
        
        if (abs(res_DimB)>1e-6 || abs(res_DimD)>1e-6)
            error('Error between the computations is unsustainable')
        end

    end


   

    if nargin==0
        if(abs(tau_D+tau_B-check_Eff)>1e-6)
            error('There is a problem with Compute Effective Stress D')
        else
            disp('Nothing to confess, father. I am without sin')
        end
    end
    
    % Effective stress
    tau_eff = tau_B+tau_D;
 if non_linear_um == 1 
        varargout{1} = tau_eff;
        varargout{2} = tau_B; 
        varargout{3} = tau_D; 
        varargout{4} = ID; 
        varargout{5} = eta_um; 
    else
        varargout{1} = tau_eff;
        varargout{2} = tau_B; 
        varargout{3} = tau_D; 
        varargout{4} = ID; 

 end
end