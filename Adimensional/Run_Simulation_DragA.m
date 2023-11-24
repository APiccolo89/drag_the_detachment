function Testdata = Run_Simulation_DragA(ID_A,nlm)
%Compute the D,t of the slab detachment.
%======================================================================
% Output Parameter:
% Testdata => structure containing all the solution output that are then
% used for the data reduction
%======================================================================
% Input Parameter:
% ID_A => Initial data structure containing all the rheological parameter
% of the slab and upper mantle adimensional
% Short description:
% => Create a function handle to introduce into ode15i (variable method to
% solve non linear implicit function)
% => use and do plot
% =>save the associated vector in a structure and compute the next solution
%======================================================================
%
Benchmark = 0 ;
if(nargin==0)
    % Benchmark = 1;
    ID_ = Compute_slab_characteristics(NaN);
    ID_A  = ID_.ID_A;
    nlm   = Problem_type.Linear;
end
success = 0; % Flag to warn if the test succeded or not
iteration = 0; % flag that warn how many attempt of doing the test has been made

% 1) Assuming that drag force is not active
dDdt0 = -0.5.*ID_A.dDdtB;
% Create function handle
Funf_wi = @(t,x,xp0) compute_dragODEA(x,xp0,ID_A,nlm,0);
% Most of the time the tests are not working because the solver
% overshoot the last timestep. In order to verify this, I introduced a
% cycle while together with a try and catch routine. Not the optimal
% solution. But it is the best that I come out atm.
max_step=1e-3;
step = 0.5*1e-2;
max_val = 1e-6;
% Set the option for resolving the system of equation
options = odeset('RelTol',1e-12,'NormControl','on','Events',@(t,x,xp0) det_EV(t,x,xp0,ID_A),'MaxStep',1e-2);
% resolve the system
while success == 0
    try % See if it works
        [t,D,te,De,ie] = ode15i(Funf_wi,[0 20],ID_A.D0,dDdt0,options);
        % If it works change the flag
        success  = 1;
        % Post process the test
        [Testdata]=postprocess_data(t,D,ID_A,te,De,ie,0,0,nlm);
        % Save the information that the test is not failed
        Testdata.Failed = 0;
        % Save the relative tollerance
        Testdata.option.Tol = 1e-12;
        % Save the max step introduced
        Testdata.option.max_step = max_step;
        % Save the number of iteration
        Testdata.option.iteration = iteration;
        % Tell me something nice
       % disp(['Test have been tried',num2str(iteration),'times'])
    catch
        % nope
        success = 0;
        % decrease the step of 5e-1
        max_step = max_step*0.5*1e-1;
        % add the information
        iteration = iteration+1;
        % Is the test failed? I decided that if the max step is 1e-5 it is
        % considered failed
        if max_step < max_val
            % Update the junk 
            Testdata.Failed = 1;
            Testdata.option.Tol = 1e-12;
            Testdata.option.max_step = max_step;
            Testdata.option.iteration = iteration;
           % disp(['Test have been tried',num2str(iteration),'times, BUT FAILED'])
            break; % And exit from the loop otherwise it will be continue till the end of the time or battery or civilisation, who knows?
        end
    end
end
end

function [res] = compute_dragODEA(D,dDdt,ID_A,nlm,catch_)
% Output:
% res => Residuum of non linear equation
%======================================================================
% Input :
% D    => Current thickness
% dDdt => Current time derivatives
% ID   => Initial Data structure
% Short description: recast the equation dDdt = -D(B(tau)^n+B(tau)) in
% terms of residuum. This function is the right form to use with ODE15.
% and cast the problem of necking in term of optimization.
%======================================================================
% Compute the effective stress

[tau_eff,tau_B,tau_D]         = Compute_Effective_StressA(D,dDdt,ID_A,nlm);

% Compute the strain rate
[epsilon_eff,eps_dif,eps_dis] = Compute_StrainA(ID_A,tau_eff);
% Necking rate equation
res                           = -D*(epsilon_eff)-dDdt;
end
% Event detection (not difficult, though, i just copied the matlab help
% page (which has been copied also by Marcel in his code)
function [position,isterminal,direction] = det_EV(t,y,yp0,ID)
position = y(1)-0.1; % The value that we want to be zero
% To be precise, each time step y(1) is the actual thickness. so,
% D(t)-0.1*D0 must be 0 to stop the simulation and having an event.
% Additional mistake that I did the function handle must incorporate the
% derivative too.
isterminal = 1;  % Halt integration
direction = 0;   % The zero can be approached from either direction
end
