function [tau_eff,tau_B,tau_D] = Compute_effective_StressD(D,dDdt,ID)
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
    end
    % Buoyancy stress computed tau_B = F_B/2/D; 
    tau_B = (ID.s0*ID.D0)./D;
    % Drag force related stress compute formulation of Bercovici et al 2015 
    tau_D = (2*ID.etaum*ID.alpha*(ID.D0^2./D.^2).*dDdt.*ID.len)./2./D;
    
    if nargin==0
        if(abs(tau_D+tau_B-check_Eff)>1e-6)
            error('There is a problem with Compute Effective Stress D')
        else
            disp('Nothing to confess, father. I am without sin')
        end
    end
    
    % Effective stress
    tau_eff = tau_B+tau_D;

end