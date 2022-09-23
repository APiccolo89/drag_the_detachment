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


function [etaum,tau_D,ID] = compute_iteration(ID,dDdt,D,eta_um_)

    tau_D_  = abs(2*(1/ID.Df_UM).*eta_um_*ID.alpha*(ID.D0^2./D.^2).*dDdt.*ID.len)./2./D;
    lim = 1e-6;
    res = 1.0;
    it = 0;
    loop = 1; 
    while loop == 1
        eta_um_0 = (ID.etaum)./(1+ID.Df_UM*(abs(tau_D_./ID.s0)^(ID.n-1)));
        tau_D_0  = abs(2*eta_um_0*ID.alpha*(ID.D0^2./D.^2).*dDdt.*ID.len)./2./D;
        res = abs(tau_D_0-tau_D_);
        if eta_um_0 < 1e13
            error('C')
        end
        tau_D_=tau_D_0;
        it = it+1; 
        if res < lim || it > 100
             loop = 0; 
        end
      
    end
    etaum = eta_um_0;
    ID.eta_um_ = etaum; 
    tau_D = tau_D_0*sign(dDdt);

    end