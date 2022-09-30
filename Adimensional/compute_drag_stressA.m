function [tau_D,Lambda] = compute_drag_stressA(ID,D,dDdt)
%Input:
%==========================================================================
% ID initial data structure
% D  the actual thickness
% dDdt the actual necking velocity
%==========================================================================
%Output:
%==========================================================================
%tau_D converged tau_D
%converged,eta_um
%
%==========================================================================
% Compute the drag force with an itera  tive algoritm (rel tol 1e-8)
% Guess assuming that the initial etaum is the diffusion viscosity
% compute tau_D
%==========================================================================
if nargin == 0
[ID_] = Compute_slab_characteristics();
ID = ID_.ID_A; 
dDdt = -1; 
D = 1.0;  
end




if length(dDdt)>1
    n_points = length(dDdt);
    Lambda = zeros(n_points,1);
    tau_D = zeros(n_points,1);
    for i=1:n_points
        [Lambda(i), tau_D(i)] =  compute_drag_stress_etaA(ID,dDdt(i),D(i));
    end
else

    [Lambda, tau_D] =  compute_drag_stress_etaA(ID,dDdt,D);

if nargin == 0 
    [tau_D_D, eta_um,DIM] = compute_drag_stress(); 
    tau_D_D_ = tau_D_D/DIM.s0;
    res_stress = tau_D_D_-tau_D; 
    % compute Lambda: 
    Lambda_D   = DIM.Lambda/(1+(DIM.Df_UM)*(abs(tau_D_D/DIM.s0))^(DIM.n-1));
    res_lambda = Lambda_D-Lambda; 
    disp('  Residuum between computed by adimensional and dimensional is: ')
    disp(['     tau_D(D)-tau_D(A) = ', num2str(res_stress,3)])
    disp(['     Lambda(D)-Lambda(A) = ', num2str(res_lambda,3)])
    disp(['===================================================================='])

end




end

end


function [Lambda,tau_D] = compute_drag_stress_etaA(ID,dDdt,D)
%=====================================================================%
% Input: ID   : structure containing the initial data (to improve)    %
%        dDdt : the converged dDdt                                    %
%        D    : actual thickness                                      %
%======================================================================
% Output: Lambda : Lambda becomes time dipendent: lambda0/(1+(t_D/t_B))
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
eps_A = abs((ID.alpha/ID.s)*(ID.D0/D)^2*dDdt);
eta_n = (2*ID.B_n_um^(1/ID.n)*eps_A^((ID.n-1)/ID.n))^(-1);
eta_eff_g =(1/ID.etaum+1/eta_n)^-1;
% Compute the new initial guess Lambda
Lambda_g = (eta_eff_g/ID.etaS)*(ID.L0*ID.alpha)/ID.s;
% fzero
tau_g    = 2*eta_eff_g*eps_A;
fun      = @(x) f_zero_L_T(x,dDdt,D,ID,eps_A);
tau_D    = fzero(fun,tau_g);
Lambda   = ID.Lambda/(1+ID.Df_UM*(abs(tau_D))^(ID.n-1));
tau_D    = tau_D*sign(dDdt);

end



function [res] = f_zero_L_T(tau_D,dDdt,D,ID,eps_A)
%====================================================================%
% Cast the problem of finding the root as optimization parameter
% use fzero and retrieve the real stress
% use the real stress to compute the new lambda.
%====================================================================%
%Lambda_r = ID.Lambda/(1+ID.Df_UM*(abs(tau_D))^(ID.n-1));
eps = ID.B_n_um*(abs(tau_D)^ID.n)*(1+(1/ID.Df_UM)*(abs(tau_D))^(1-ID.n));
res = eps_A-eps;
end