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
%Picard iteration scheme initial guess -> compute a value, use the value as
%initial guess -> compute a new one till abs(x(n-1)-xn)<1e-10
%==========================================================================
% Compute the drag force with an iterative algoritm (rel tol 1e-8)
% Guess assuming that the initial etaum is the diffusion viscosity
% compute tau_D
%==========================================================================
if length(dDdt)>1
    n_points = length(dDdt);
    Lambda = zeros(n_points,1);
    tau_D = zeros(n_points,1);
    for i=1:n_points
    [Lambda(i), tau_D(i)] =  compute_picard(ID,dDdt(i),D(i)); 
    end
else

 [Lambda, tau_D] =  compute_picard(ID,dDdt,D); 


end

end


    function [Lambda,tau_D_] = compute_picard(ID,dDdt,D)

    tau_B = (ID.D0./D);
    Lambda_0 = ID.Lambda./(1+ID.Df_UM*(1.0)^(ID.n-1));
    tau_D_ = abs(ID.Lambda.*tau_B.^3.*dDdt);
    lim = 1e-6;
    res = 1.0;
    it = 0;
    loop = 1; 
    while loop == 1
        Lambda_0 = ID.Lambda./(1+ID.Df_UM*(tau_D_)^(ID.n-1));
        tau_D_0  = abs(Lambda_0.*tau_B.^3.*dDdt);
        res = abs(tau_D_0-tau_D_);
        tau_D_=tau_D_0;
        it = it+1; 
        if res < lim || it > 100
             loop = 0; 
        end
      
    end
    Lambda = Lambda_0;
    tau_D_ = tau_D_0*sign(dDdt);

    end
    