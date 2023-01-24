function [UPPER_MANTLE,SLAB] = main_function_Real(t0, T_Mean, Tp,Pr, D0,L0,Vnv,Vdv,UM,S)

% compute mantle averages along L 

% compute slab rheology as a function average T
[MMVn,MMVd,MML0]  = meshgrid(Vnv,Vdv,L0);
eta0DMP = (MMVn./MMVn).*0.0;
eta0MP  = (MMVn./MMVn).*0.0;
xiumP   = (MMVn./MMVn).*0.0;

for i=1:length(MMVn(:))
    l0 = MML0(i);
    Vn = MMVn(i);
    Vd = MMVd(i);
    [eta0DMP(i), eta0MP(i), xiumP(i)] = compute_effective_viscosities_Mantle(t0,l0,Vn,Vd,Tp,Pr,UM);
end
% Store in structure Upper mantle
UPPER_MANTLE.eta0DMP = eta0DMP;
UPPER_MANTLE.eta0MP = eta0MP;
UPPER_MANTLE.xiumP  = xiumP;
UPPER_MANTLE.Vd     = MMVd; 
UPPER_MANTLE.Vn     = MMVn;
UPPER_MANTLE.L0     = MML0; 

[MMvnS,MMvdS,MMT_av] = meshgrid(Vnv,Vdv,T_Mean);
eta0DS = (MMvnS./MMvnS).*0.0;
eta0S  = (MMvnS./MMvnS).*0.0;
xiuS   = (MMvnS./MMvnS).*0.0;

for i = 1:length(MMvnS(:))
T_avg = MMT_av(i); 
VnS    = MMvnS(i); 
VdS    = MMvdS(i); 
[eta0S(i),eta0DS(i)]=compute_effective_viscosity_slab(VnS,VdS,t0,Pr,T_avg,S);
xiuS(i) = eta0DS(i)./eta0S(i);
end
% Store Slab structure 

SLAB.eta0DS = eta0DS;
SLAB.eta0S = eta0S;
SLAB.xiuS  = xiuS;
SLAB.Vd     = MMvdS; 
SLAB.Vn     = MMvnS;
SLAB.MMT_av     = MMT_av; 

end


function [eta0MDP,eta0MP,xium_P] = compute_effective_viscosities_Mantle(t0,L0,Vn,Vd,Tp,Pr,UM)  
 
% Constants
drho = (t0*2)/(9.81*L0);
[Cd,Cn,phi] = UM.Compute_Cd_Cn(Pr,Vn,Vd,Tp); 
% Compute Reference viscosity 
eta0M = compute_reference_viscosity(UM.Bn,UM.En,Vn,UM.n,Tp,Pr,t0,UM.R);
eta0DM = compute_reference_viscosity(UM.Bd,UM.Ed,Vd,1,Tp,Pr,t0,UM.R);
xium_0 = eta0DM/eta0M; 
% Compute integral Upper mantle 
Exn = @(x) compute_exponential(Cn,x,phi,3300*UM.g);
Exd = @(x) compute_exponential(Cd,x,phi,3300*UM.g);
intexn = (integral(Exn,0,L0))./L0;
intexd = (integral(Exd,0,L0))./L0;
eta0MP = eta0M*intexn;
eta0MDP = eta0DM*intexd; 
xium_P = eta0MDP/eta0MP; 
if isnan(xium_P)
    bla = 0; 
end
disp('=====================================================================')
disp(['Average reference viscosity UM integrated along the slab = '])
disp(['log10(eta0) = ',num2str(log10(eta0MP),4)])
disp(['log10(eta0D) = ', num2str(log10(eta0MDP),4)])
disp(['xiUM0 = ' num2str(log10(xium_0)), ' and xiUMP =', num2str(log10(xium_P))])
disp('=====================================================================')

end

function  [eta0SP,eta0DSP] = compute_effective_viscosity_slab(Vn,Vd,t0,Pr,T_Mean,S)
[Cd,Cn,phi] = S.Compute_Cd_Cn(Pr,Vn,Vd,T_Mean); 
% Compute Reference viscosity slab 
Tavg =  T_Mean;%compute_average_T_slab(D0,age,Tp);
if Tavg > 800
bla = 0 
end
eta0SP = compute_reference_viscosity(S.Bn,S.En,Vn,S.n,Tavg,Pr,t0,S.R);
eta0DSP = compute_reference_viscosity(S.Bd,S.Ed,Vd,1,Tavg,Pr,t0,S.R);

end



function [exp_] = compute_exponential(C,l,phi,w)
    exp_= exp(-C.*((l.*w)./(1+phi.*l*w)));
end
