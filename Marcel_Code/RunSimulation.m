% run a 0D necking simulation
%------------------------------
function [t,y,te,ye,ie] = RunSimulation(D,S,r_sdis0,deltaQg,deltaQdif,n,m,q,T0,r_init,d_neck)

if nargin == 0
    close all
    % parameter setting
    %--------------------------
    D           = 0; % D==0 switches off grain size evolution
    S           = 0; % S==0 switches off thermal evolution
    r_sdis0     = 0.0; %? stress
    deltaQg     = 0.56; % Thermal parameter ?
    deltaQdif   = 0.56; % Thermal parameter ?
    n           = 3.5;  % stress exponent
    m           = 0;    % grain size exponent
    q           = 0;    % grain exponent growth
    
    d_neck = 1e-4;    %? final neck size
    
    % initial parameters
    %--------------------------
    T0 = 57;%5.3*10^5/8.3144 / 1300; % this is the nondimensional parameter!!! not the initial nondimensional temperature perturbation!!!
    r_init = 0; % determines whether we start in diffusion or dislocation creep,r_init>1--> dislocation
    
end

T_init = 0; % by definition
d_init = 1; % by definition

% add some scaling to the time, as things are sometimes a bit too fast for the solver
tscale = 1; % this is completely arbitrary, check if it really works

% initial solution
Sol0        = [r_init;T_init;d_init];

% time span to simulate
tspan       = [0 1].*tscale; % should be detached at 1 for Newtonian no damage case
% function handle
ODEfun      = @(t,x) Compute_dSoldt(t,x,D,S,r_sdis0,deltaQg,deltaQdif,n,m,q,T0,tscale);
%%% I've never used a function handle before, so better to understand
%%% clearly: He passes a function with two input parameters to a function
%%% that is composed by three subfunction, then, it passes to the solver
% necking event

% options
options = odeset('RelTol',1e-8,'NonNegative',1:3,'NormControl','on', 'Events',@(t,x) AllEvents(t,x,d_neck,deltaQdif,n,m,T0),'InitialStep',1e-10,'Refine',4);

% run
% ode23tb
% ode15s
[t,y,te,ye,ie] = ode23tb(ODEfun, tspan, Sol0,options);


te = te./tscale;
t = t./tscale;

% process events
ind_end = find(ie==1);

if ~isempty(ind_end)
    t_neck = te(ind_end);
end

iloop = 1;
while isempty(ind_end) % this should happen if there is significant necking
    if iloop>5
       %tscale = tscale*100; 
       
       if tscale >1e14
          bla = 1; 
       end
    end
    
    [t2,y2,te2,ye2,ie2] = ode23tb(ODEfun, tspan, y(end,:),options);
    
    % rescale
    t2   = t2./tscale;
    te2 = te2./tscale;
    
    % add events
    te = [te;te2+t(end)];
    ye = [ye;ye2];
    ie = [ie;ie2];
    
    ind_end = find(ie==1);
    if ~isempty(ind_end)
        % necking time
        t_neck = (te(ind_end) + t(end));
    end
    % piece the solution together
    t = [t;t2(2:end)+t(end)]; %remove the scale
    y = [y;y2(2:end,:)]; 
    
    iloop = iloop+1;
end

if nargin == 0
    %% postprocessing:
    % separate solution components
    r = y(:,1);
    T = y(:,2);
    d = y(:,3);
    
    % - compute strain rates of each mechanisms
    % - compute thermal and grain size contribution to strain rate
    [edis,edis_T,edif,edif_T,edif_r] = ComputeStrainRates(r,T,d,deltaQdif,n,m,T0);
    
    ind_disdif = find(ie==2);
    ind_difdis = find(ie==3);
    ind_Tr     = find(ie==4);
    ind_rT     = find(ie==5);
    
    %% plot
    figure(1)
    subplot(231) % slab thickness
    hold on
    plot(t,d)
%     p(1) = plot(te(ind_disdif),ye((ind_disdif),3),'ro');
%     p(2) = plot(te(ind_difdis),ye((ind_difdis),3),'go');
%     p(3) = plot(te(ind_Tr),ye((ind_Tr),3),'ks');
%     p(4) = plot(te(ind_rT),ye((ind_rT),3),'bv');
%    legend(p,{'dis2dif';'dif2dis';'T2r';'r2T'},'Location','southwest')
    set(gca,'YScale','log')
    
    subplot(232) % temperature perturbation
    hold on
    plot(t,T)
%     p(1) = plot(te(ind_disdif),ye((ind_disdif),2),'ro');
%     p(2) = plot(te(ind_difdis),ye((ind_difdis),2),'go');
%     p(3) = plot(te(ind_Tr),ye((ind_Tr),2),'ks');
%     p(4) = plot(te(ind_rT),ye((ind_rT),2),'bv');
%    legend(p,{'dis2dif';'dif2dis';'T2r';'r2T'},'Location','southwest')
    set(gca,'YScale','log')
    subplot(233) % roughness
    hold on
    plot(t,r)
%     p(1) = plot(te(ind_disdif),ye((ind_disdif),1),'ro');
%     p(2) = plot(te(ind_difdis),ye((ind_difdis),1),'go');
%     p(3) = plot(te(ind_Tr),ye((ind_Tr),1),'ks');
%     p(4) = plot(te(ind_rT),ye((ind_rT),1),'bv');
%    legend(p,{'dis2dif';'dif2dis';'T2r';'r2T'},'Location','southwest')
    set(gca,'YScale','log')
    
    subplot(234)
    hold on
    plot(t,edis)
    plot(t,edif)
    legend({'dis';'dif'})
    set(gca,'XScale','linear','YScale','log')
    subplot(235)
    hold on
    plot(t,edis./edif)
    set(gca,'XScale','linear','YScale','log')
    subplot(236)
    hold on
    plot(t,edif_T./edif_r)
%     p(1) = plot(te(ind_Tr),ye((ind_Tr),3),'ro');
%     p(2) = plot(te(ind_difdis),ye((ind_difdis),3),'go');
%     p(3) = plot(te(ind_Tr),ye((ind_Tr),3),'ks');
%     p(4) = plot(te(ind_rT),ye((ind_rT),3),'go');
%    legend(p,{'T2r';'r2T'},'Location','southwest')
    set(gca,'YScale','log')
    bla = 1;
end



%% EVENT FUNCTION FOR NECKING

function [position,isterminal,direction] = NeckEvent(t,y,d_neck)
position = y(3)-d_neck; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction

function [position,isterminal,direction] = Dis2DifEvent(t,y,d_neck,deltaQdif,n,m,T0)
[edis,edis_T,edif,edif_T,edif_r] = ComputeStrainRates(y(1),y(2),y(3),deltaQdif,n,m,T0);
position = edis/edif-1; % The value that we want to be zero
isterminal = 0;  % Halt integration 
direction = 1;   % The zero can be approached from the positive side

function [position,isterminal,direction] = Dif2DisEvent(t,y,d_neck)
[edis,edis_T,edif,edif_T,edif_r] = ComputeStrainRates(y(1),y(2),y(3),deltaQdif,n,m,T0);
position = edis/edif-1; % The value that we want to be zero
isterminal = 0;  % Halt integration 
direction = -1;   % The zero can be approached from either direction

function [position,isterminal,direction] = AllEvents(t,y,d_neck,deltaQdif,n,m,T0)
[edis,edis_T,edif,edif_T,edif_r] = ComputeStrainRates(y(1),y(2),y(3),deltaQdif,n,m,T0);

position = [y(3)-d_neck; edis/edif-1; edis/edif-1;edif_T./edif_r-1;edif_T./edif_r-1];
isterminal = [1; 0; 0;0;0];
direction = [0; -1; 1;-1;1];


