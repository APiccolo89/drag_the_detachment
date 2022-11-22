function Testdata = Run_Simulation_Drag(ID)
    %Compute the D,t of the slab detachment. 
    % Output Parameter: 
    % Testdata => structure containing all the solution output that are then
    % used for the data reduction 
    % Input Parameter: 
    % ID => Initial data structure containing all the rheological parameter
    % of the slab and upper mantle
    % Short description:
    % => Create a function handle to introduce into ode15i (variable method to
    % solve non linear implicit function)
    % => use and do plot
    % =>save the associated vector in a structure and compute the next solution
    %%
    %Main Function
    
    % Compute the initial guess of the dDdt 
    % 1) Assuming that drag force is not active
    dDdt0 = compute_guess_ddt(ID);
    % Create function handle 
    Funf_wi = @(t,x,xp0) compute_dragODE(x,xp0,ID);
    % Set the option for resolving the system of equation
    options = odeset('RelTol',1e-8,'NormControl','on');
    % resolve the system
    [t,D] = ode15i(Funf_wi,[0 ID.tc],ID.D0,dDdt0,options);
    % Normalize the thickness vector
    D_norm = D/ID.D0;
    % save relevant data of the simulation:
    Testdata(1,:)=t/ID.tc/n; %time vector divided by the detachment timescale 
    Testdata(2,:)=D_norm;
end

                string_ID = {'B_d','B_n','s0','n','eta0','Df','drho','D0','l0','etaum'};

function [dDdt0] = compute_guess_ddt(ID)
    % Output:
    % [dDdt0] => initial buoyancy related necking velocity
    % Input: 
    % ID => Initial data structure
    % 
    % Compute the initial dDdt assuming that no drag forces are active
    
    % Buoyancy stress computed using F_B formulation
    tau = (ID.drho*ID.D0*ID.l0*9.81)/2/D;
    
    % Initial necking velocity
    dDdt0 = (-D*(B_n*tau^n+B_d*tau));
end


function [res] = compute_dragODE(D,dDdt,ID)
    % Output: 
    % res => Residuum of non linear equation 
    % Input : 
    % D    => Current thickness
    % dDdt => Current time derivatives
    % ID   => Initial Data structure 
    % Short description: recast the equation dDdt = -D(B(tau)^n+B(tau)) in
    % terms of residuum. This function is the right form to use with ODE15.
    % and cast the problem of necking in term of optimization. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the effective stress 
    [tau_eff] = compute_the_effective_stress(D,dDdt,ID,0);
    res = -D*(ID.B_n*(tau_eff)^n+ID.B_d*(tau_eff))-dDdt;
end


function [tau_eff] = compute_the_effective_stress(D,dDdt,ID)
    % Output: 
    % tau_eff  = effective stress 
    % Input :
    % D  => Actual thickness of the slab
    % dDdt => Actual rate of necking
    % ID => Initial data structure 
    %
    % Short description: 
    % Compute the stress related to the buoyancy (tau_B) 
    % Compute the stress associated with the drag force acting in the lateral 
    % boundary of the slab (tau_D)
    % tau_eff = tau_B+tau_D
    l = ID.l0/1000e3; % initial length divded by the size of the model (or initial thickness of the slab)
    % Buoyancy stress computed tau_B = F_B/2/D; 
    tau_B = (ID.drho*ID.D0*ID.l0*9.81)/2/D;
    % Drag force related stress compute formulation of Bercovici et al 2015 
    tau_D = (2*ID.etaum*alpha*(ID.D0^2/D^2)*dDdt*l)/2/D;
    % Effective stress
    tau_eff = tau_B+tau_D;
end


