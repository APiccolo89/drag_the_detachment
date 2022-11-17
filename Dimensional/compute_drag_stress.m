function [tau_M,etaum,ID] = compute_drag_stress(ID,D,dDdt)
%Input:
% UPDATE: I'm Idiot and I was confusiong the drag stress that i compute
% here for the one that actually is computed for the slab. 
% I.e. IN THIS FUNCTION I ESTIMATE THE AVERAGE VISCOSITY WITHIN A COUETTE
% LIKE FLOW WHOSE DIMENSION IS THE LENGHT OF THE SLAB!!!!! AND I NEED TO DO
% MAJOR UPDATE FOR THE CODE AND JUST TRACK THE CHANGE!!!!
%
%==========================================================================
% ID initial data structure
% D  the actual thickness
% dDdt the actual necking velocity
%==========================================================================
%Output:
%==========================================================================
%tau_M converged tau_M converged mantle stresses assuming that the average
%stress within the couette flow arises from the average strain within the
%couette flow
%converged,eta_um
% iteration scheme initial guess -> compute a value, use the value as
%initial guess -> compute a new one till abs(x(n-1)-xn)<1e-10
%==========================================================================
% Compute the drag force with an iterative algoritm (rel tol 1e-8)
% Guess assuming that the initial etaum is the diffusion viscosity
% compute tau_D
% Fun fact: I need to keep track of the viscosity of the upper mantle, as
% the dDdt is evolving as a consequence of the viscosity. 
%
%
%==========================================================================
if nargin == 0
    [ID_] = Compute_slab_characteristics();
    ID = ID_;  
    D = ID.D0;
    dDdt = -D*(ID.ec);
end



if length(dDdt)>1
    n_points = length(dDdt);
    etaum = zeros(n_points,1);
    tau_M = zeros(n_points,1);
    for i=1:n_points
    [etaum(i), tau_M(i)] =  compute_drag_stress_etaD(ID,dDdt(i),D(i)); 
    end
else

 [etaum, tau_M] =  compute_drag_stress_etaD(ID,dDdt,D); 

 if nargin == 0 
    [tau_M_,Lambda_] = compute_drag_stressA(ID.ID_A,D/ID.D0,dDdt/(D*(ID.ec))); 
    res_stress = tau_M_-tau_M/ID.s0; 
    % compute Lambda: 
    Lambda_A   = ID.Lambda/(1+(ID.Df_UM)*(abs(tau_M_))^(ID.n-1));
    Lambda_D   = (etaum/((1/ID.eta0NS+1/(ID.eta0DS))^(-1)))*((ID.l0*ID.alpha)/(2*ID.s));
    res_lambda = Lambda_A-Lambda_D; 
    disp(['===================================================================='])
    disp('  Residuum between computed by dimensional and adimensional is: ')
    disp(['     tau_D(A)-tau_D(D) = ', num2str(res_stress,3)])
    disp(['     Lambda(A)-Lambda(D) = ', num2str(res_lambda,3)])
    disp(['===================================================================='])

end

end

end


function [eta_um,tau_M] = compute_drag_stress_etaD(ID,dDdt,D)
%=====================================================================%
% Input: ID   : structure containing the initial data (to improve)    %
%        dDdt : the converged dDdt                                    %
%        D    : actual thickness                                      %
%======================================================================
% Output: eta_um : viscosity of the upper mantle                      %
%         tau_d_ : the drag stress.                                   %
%======================================================================
% Alg : a) Compute guess stress using dDdt.                           %
%       b) Compute new stress and then viscosity with fzero           %
%       c) Output viscosity and tau_D                                 %
%======================================================================
%Lengthy explanation: Initially I tried to compute the guess using the
%stress. In reality I was complicating my life. I could have used the
%effective strain rate (vz/s) provided by the formulation of Bercovici.
%I compute the dislocation creep guess using the initial strain rate.
%then I compute the effective viscosity and compute a lambda guess.
%After that, I compute Lambda using a stress based definition, which is
%more safe (eps_tot = eps_dif+eps_dis) while tau is the same for all
%the mechanism and to be consistent with my derivation.
%======================================================================
% Prepare the initial value
% compute the initial strain rate

eps_A = abs(ID.alpha*(ID.D0^2./D.^2).*dDdt.*1/ID.l0);
eta_n = (2*ID.B_n_um^(1/ID.n)*eps_A^((ID.n-1)/ID.n))^(-1);
eta_eff_g =(1/ID.eta0DM+1/eta_n)^-1;
%options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
tau_g    = (2*eta_eff_g*eps_A);                               %Here there is a mistake, perhaps correct intuition, but i need a way to find out 
fun      = @(x) f_zero_L_TD(x,dDdt,D,ID,eps_A);
tau_M    = fzero(fun,abs(tau_g));
eta_um   = ID.eta0DM/(1+ID.Df_UM*(abs(tau_M)/ID.s0)^(ID.n-1));
tau_M    = tau_M*sign(dDdt); 
end



function [res] = f_zero_L_TD(tau_D,dDdt,D,ID,eps_A)
%====================================================================%
% Cast the problem of finding the root as optimization parameter
% use fzero and retrieve the real stress
% use the real stress to compute the new etaum.
%====================================================================%
eps = ID.B_n_um*(abs(tau_D)^ID.n)*(1+(1/ID.Df_UM)*(abs(tau_D)/ID.s0)^(1-ID.n));
%eta_n = (2*ID.B_n_um*abs(tau_D)^((ID.n-1)))^(-1);
%eta_eff_g =(1/ID.etaum+1/eta_n)^-1;
res = eps-eps_A;   
% (eps_A-(1/(2*ID.etaum)*abs(tau_D)+ID.B_n_um*abs(tau_D)^(ID.n)));
end