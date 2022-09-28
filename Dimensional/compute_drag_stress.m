function [tau_D,etaum,ID] = compute_drag_stress(ID,D,dDdt)
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
if length(dDdt)>1
    n_points = length(dDdt);
    etaum = zeros(n_points,1);
    tau_D = zeros(n_points,1);
    for i=1:n_points
    if i == 1 
        eta_um_ = ID.etaum; 
    else
        eta_um_ = etaum(i-1); 
    end
    [etaum(i), tau_D(i),ID] =  compute_iteration(ID,dDdt(i),D(i),eta_um_); 
    end
else

 [etaum, tau_D,ID] =  compute_iteration(ID,dDdt,D,ID.eta_um_); 


end

end


function [eta_um,tau_D,ID] = compute_iteration(ID,dDdt,D,eta_um_)
    s = ID.s;               % size of the convection cell 
    alpha = ID.alpha;       % alpha 
    B_n_um  =   (ID.s0^(1-ID.n))/ID.etaum;
    epsilon_tot = abs(alpha*(ID.D0^2./D.^2).*dDdt.*1/s);
    eta_n_ = 0.5*(1/(ID.B_n_um^(1/ID.n)*epsilon_tot^((ID.n-1)/ID.n))); 
    eta_d_ = ID.etaum; 
    eta_um = (1/eta_n_+1/eta_d_)^(-1);
    tau_D  = abs(2*eta_um*ID.alpha*(ID.D0^2./D.^2).*dDdt.*ID.len)./2./D;
    end