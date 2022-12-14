function Testdata = Run_Simulation_Drag(ID,Benchmark,nlm)
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
    if(nargin==0)
        addpath ../Adimensional/
        Benchmark = 1;  
        nlm = Problem_type.Linear;
        ID = Compute_slab_characteristics(NaN);
    end
    
    % Compute the initial guess of the dDdt 
    % 1) Assuming that drag force is not active
    dDdt0 = -(ID.tc)^(-1)*ID.D0;
    % Create function handle 
    Funf_wi = @(t,x,xp0) compute_dragODE(x,xp0,ID,1,nlm);
    % Set the option for resolving the system of equation
    options = odeset('RelTol',1e-12,'NormControl','on','Events',@(t,x,xp0) det_EV(t,x,xp0,ID));
    % resolve the system
    [t,D,te,De,ie] = ode15i(Funf_wi,[0 20*ID.tc],ID.D0,dDdt0,options);

    
    

    % Normalize the thickness vector
    % Function to post process the stress, strain and so forth
    % place holder
    % save relevant data of the simulation:
    [Testdata]=postprocess_data(t,D,ID,te,De,ie,1,Benchmark,nlm); 
    
    if nargin == 0 
    % run the twin test adimensional and check the residuum between the two
    % functions. {Is giving allmost the same results?)
    Testdata_a = Run_Simulation_DragA(ID.ID_A);
    % Interpolate the results and plot them: 
    [Intp_data]=Interpolation_routinesAD(Testdata,Testdata_a,Benchmark);
    plot_benchmark(Intp_data,ID,Testdata.D_norm,1.0)
    end
end

function [res] = compute_dragODE(D,dDdt,ID,Benchmark,nlm)
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
    [tau_eff,tau_B,tau_D,ID] = Compute_effective_StressD(D,dDdt,ID,Benchmark,nlm);
    [eps_eff,eps_dif,eps_dis] = Compute_StrainD(ID,tau_eff,Benchmark);
    D_ = D/ID.D0;
    eps_ = eps_eff/ID.ec;
    dDdt_ = dDdt/(ID.D0/ID.tc);
    res = -D_*(eps_)-dDdt_;
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
