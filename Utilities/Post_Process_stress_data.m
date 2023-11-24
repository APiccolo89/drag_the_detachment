function [tau_m,eta_mantle,xdisl_art] = Post_Process_stress_data(Bn,Xi,eta,theta,n)
Gamma  = (1-n)./n;

for i = 1:length(Xi(:))
    A = cputime;
    % Find the stress using fzero    
    % Since I am effectively mentally retarded, i was using my usual
    % approach for finding the dislocation creep rheology, neglecting the
    % small and pesky details that the characteristic strain rate is
    % defined using a reference stress, so, the strain rate of the upper
    % mantle is not the same of the one that you might find in the slab. 
    % Computing the bn 
    %Bn = 1./(2.*eta(i)./Xi(i));
    % Computing bd
    %Bd = 1./(2.*eta(i));
    % Computing the strain rate indipendent part
    eta0=0.5 * Bn(i)^(-1/n);
    % Computing the dislocation creep
    eta_disl = eta0.*theta.^Gamma;
    % Computing the maximum biscosity
    eta_max = min(eta(i),eta_disl);
    %Computing the minimum viscosity
    eta_min = (1./eta_disl+1./eta(i)).^-1; 
    %Computing the mean viscosity
    eta_mean = 0.5.*eta_min+0.5.*eta_max;
    %Computing the guess of the stress
    tau_guess=[(eta_min.*2.*theta)+(eta_max.*2.*theta)]./2;
    %Computing the stress
    tau_m(i) = fzero(@(x)(x/(2 * eta(i))) + ((x/(2*eta0))^n) - theta, tau_guess);
    % Computing the real viscosity
    eps_disl = Bn(i).*tau_m(i).^n;
    %Computing the partition
    xdisl_art(i) = eps_disl./theta;
    %Computing the real effective viscosity
    eta_mantle(i) = 2.*theta.*tau_m(i);
    B = cputime;
   
    disp(['time =',num2str(B-A), 's']);
end

end

