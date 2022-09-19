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
    dDdt0 = -(ID.tc)^(-1)*ID.D0;
    % Create function handle 
    Funf_wi = @(t,x,xp0) compute_dragODE(x,xp0,ID);
    % Set the option for resolving the system of equation
    options = odeset('RelTol',1e-8,'NormControl','on','Events',@(t,x,xp0) det_EV(t,x,xp0,ID));
    % resolve the system
    [t,D,te,De,ie] = ode15i(Funf_wi,[0 1000*ID.tc],ID.D0,dDdt0,options);
    % Normalize the thickness vector
    % Function to post process the stress, strain and so forth
    % place holder
    % save relevant data of the simulation:
    D_cen = 0.5*(D(1:1:end-1)+D(2:1:end));
    t_c   = 0.5*(t(1:1:end-1)+t(2:1:end));
    dDdt_c = (D(2:1:end)-D(1:end-1))./(t(2:1:end)-t(1:1:end-1));
    t_B    = (ID.drho*ID.D0*ID.l0*9.81)./2./D_cen;
    t_D    = (2*ID.etaum*5.0.*(ID.D0^2./D_cen.^2).*dDdt_c*ID.l0/1e6)./2./D_cen;
  

    Testdata.time   = t_c/ID.tc; %time vector divided by the detachment timescale 
    Testdata.D_norm = D_cen/ID.D0;
    Testdata.t_det   = te/ID.tc;
    Testdata.tau(1,:) = t_B./ID.s0;
    Testdata.tau(2,:) = t_D./ID.s0;
    Testdata.tau(3,:) = (t_B+t_D)./ID.s0;
    % Find the max tau eff
    t_t_max = max(Testdata.tau(3,:));
    ind     = find(Testdata.tau(3,:)==t_t_max);
    time_t_M = t_c(ind);
    Testdata.t_t_max = t_t_max;
    Testdata.time_t_M = time_t_M./ID.tc; 
    Testdata.t_t_det  = Testdata.tau(3,end);
    %disp(['time of detachment is t/tc  ',num2str(te/ID.tc,4),' which is td_O/td_P ', num2str((ID.n*te)/ID.tc,4), '  w.r.t. analytical prediction'])
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
    [tau_eff,tau_B,tau_D] = Compute_effective_StressD(D,dDdt,ID);
    [eps_eff,eps_dif,eps_dis] = Compute_StrainD(ID,tau_eff);
    res = -D*(eps_eff)-dDdt;
end




% Event detection (not difficult, though, i just copied the matlab help
% page (which has been copied also by Marcel in his code) 
function [position,isterminal,direction] = det_EV(t,y,yp0,ID)
    position = y(1)-0.1*ID.D0; % The value that we want to be zero
    % To be precise, each time step y(1) is the actual thickness. so,
    % D(t)-0.1*D0 must be 0 to stop the simulation and having an event. 
    % Additional mistake that I did the function handle must incorporate the
    % derivative too. 
    isterminal = 1;  % Halt integration 
    direction = 0;   % The zero can be approached from either direction
end