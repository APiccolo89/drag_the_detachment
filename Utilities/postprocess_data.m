function [Testdata]=postprocess_data(t,D,ID,te,De,ie,Dim,Benchmark,nlm)
% =====================================================================
% post process the data using the same function used to compute the raw
% results.
%======================================================================
% Input:
% t   : time computed by the ODE i5
% D   : thickness evolution with time
% ID  : Initial Data Structure
% Dim : 1 means dimensional, 0, means a-dimensional
%======================================================================
% Output:
% Testdata : Data structure with the post process data.
%======================================================================
% Alg: Compute dDdt using a simple diff (diff(D)/diff(t))
%      Compute the central point of D, t
%      Post process the stress
%      Post process the strain
%      Integrate the non dimensional data in the data structure
%      Compute the detachment structure
%======================================================================
dDdt = (D(2:1:end)-D(1:end-1))./(t(2:1:end)-t(1:1:end-1));
D = 0.5*(D(1:1:end-1)+D(2:1:end));
t   = 0.5*(t(1:1:end-1)+t(2:1:end));

if nlm.islinear == 0
[t_eff,t_B,t_D,Lambda,eta_um] = Compute_Effective_StressA(D,dDdt,ID,nlm);
else
[t_eff,t_B,t_D,Lambda] = Compute_Effective_StressA(D,dDdt,ID,nlm);
eta_um = [];

end
[eps_eff,eps_d,eps_n] = Compute_StrainA(ID,t_eff);
tau = [t_B,t_D,t_eff];
eps = [eps_eff,eps_d,eps_n];
Testdata.time   = t; %time vector divided by the detachment timescale
Testdata.D_norm = D;
if isempty(te) || ~isreal(te)
    Testdata.t_det = nan; 
else
    Testdata.t_det   = te;
end
Testdata.tau(1,:) = tau(:,1);
Testdata.tau(2,:) = tau(:,2);
Testdata.tau(3,:) = tau(:,3);
Testdata.eps(1,:) = eps(:,1);
Testdata.eps(2,:) = eps(:,2);
Testdata.eps(3,:) = eps(:,3);
% Find the max tau eff
t_t_max = max(Testdata.tau(3,:));
time_t_M = t(Testdata.tau(3,:)==t_t_max);
Testdata.t_t_max = t_t_max;

Testdata.time_t_M = time_t_M;

Testdata.t_t_det  = Testdata.tau(3,end);
if nlm.islinear == 0
    Testdata.Lambda = Lambda;
    Testdata.eta_um  = eta_um;
else
    Testdata.Lambda = [];
    Testdata.etaum  = [];
end
end