%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter to test: 
% To do: 
% Generate a data structure and labelling. 

l0_v    = (300e3:50e3:600e3);                 % initial length
D0      = 80e3;                  % thickness

% Slab Rheology
eta0_v = [1e22,5e22,1e23,5e23];                   % [Pas] refernce power law viscosity slab 
Df_v   = [10];                          % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress 
n_v     =[2.0, 3.5,7.0,9.0];                       % [n.d.] pre-exponential factor
s0_v  = [100e6:50e6:200e6];                   % [Pa]  reference buoyancy stress

% Upper Mantle
etaum_v =10.^(18:0.1:21);    % [Pa.s]vector of the mantle viscosity

name_tests_v = {'D0_80km_n_2dot0_Df_EXT','D0_80km_n_3dot5_Df_EXT','D0_80km_n_7dot0_Df_EXT','D0_80km_n_9dot0_Df_EXT'};
ptsave = 'Tests_result';
if not(isfolder(ptsave))
    mkdir(ptsave)
end

for i=1:4
    name_tests = name_tests_v{i};
    n = n_v(i);

    % Function to run the ensamble of test 
    [Tests] = Run_Simulations(D0,l0_v,eta0_v,Df_v,n,s0_v,etaum_v);

    plot_results(Tests,name_tests,ptsave)
end
% Function to plot the results 



%========================================================================%
% Function utilities                                                     %
%========================================================================%
function [Tests] = Run_Simulations(D0,l0_v,eta0_v,Df_v,n,s0_v,etaum_v)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Output:
     % 
     % [Testdata] => structure containing the results of each of the
     %               simulation performed 
     % 
     % Input: 
     % 
     % D0,l0_v => Initial thickness and length (l0_v vector of data)
     %
     % eta0_v => reference power law viscosity at reference stress (tau0)
     %       n power law exponent [vector]
     % 
     % Df    =>viscosity contrast between diffusion/dislocation creep @
     %       reference stress (tau0)
     %
     % etaum_v => Viscosity of the upper mantle [vector]. 
     %
     % Alg:
     % -> loop over the data and generate combination of simulation 
     % -> construct the characteristic values 
     % -> run simulation per each combination of parameter 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     ium   = length(etaum_v);    % number of upper mantle viscosities
     iut0  = length(s0_v) ;      % number of reference buoyancy stresses
     iuet0 = length(eta0_v);     % number of reference viscosities of the slab
     iuL0  = length(l0_v)   ;    % number of tested initial length of the slab
     iuDf  = length(Df_v)   ;    % number of tested viscosity contrast at reference
     if length(n)>1 || length(D0)>1 
        error('You cannot test more than one power law exponent or initial thickness at once')
 
     end

  %loop to construct the initial data to run a simulation (to change into a
  %more performant function)
  tic
  for i=1:ium
      etaum = etaum_v(i);
      for j = 1:iut0
          s0 = s0_v(j);
        for k = 1:iuet0
            eta0 = eta0_v(k);
           for z = 1:iuL0
                l0 = l0_v(z);
                for d = 1:iuDf
                    Df = Df_v(d);
                    % Compute the characteristic for the simulation 
                    [drho,B_d,B_n,tc,ec,Psi,Lambda] = Compute_slab_characteristics(eta0,Df,n,l0,s0,D0,etaum);
                    % Create a labeling for the simulation
                    str0   = strcat('Sim_D0_',num2str(int16(D0/1e3)));
                    str1   = strcat('em_',num2str(i));
                    str2   = strcat('s0_',num2str(j));
                    str3   = strcat('eta0_',num2str(k));
                    str4   = strcat('L0_',num2str(z));
                    str5   = strcat('Df_',num2str(d));

                    l_simulation = strcat(str0,str1,str2,str3,str4,str5);
                    disp(['Simulation ',l_simulation, 'is starting'])
                    % Save into a simple data structure the initial data of
                    % the simulation
                    string_ID = {'B_d','B_n','s0','n','eta0','Df','drho','D0','l0','etaum','tc','ec','Psi','Lambda'};
                    
                    for is = 1:numel(string_ID)
                        ID.(string_ID{is}) = eval(string_ID{is});
                    end
                   
    
                    % Run a simulation with a specific combination of parameter               
                    Testdata = Run_Simulation_Drag(ID);
                    % 
                    Tests.(l_simulation) = Testdata;
                    Tests.(l_simulation).initial_data = ID; 
                    ID = [];
                    Testdata = []; 
                end
           end
        end
      end
  end
toc
end







function [drho,B_d,B_n,tc,ec,Psi,Lambda] = Compute_slab_characteristics(eta0,Df,n,l0,s0,D0,etaum)
        % Output: 
        % drho = density contrast
        % B_d  = diffusion creep pre-exponential factor
        % B_n  = dislocation creep pre-exponential factor
        % tc   = characteristic time (s)
        % ec   = characterstistc strain rate (1/s)
        % Psi  = ratio between mantle effetive viscosity and slab effective
        % viscosity
        % Lambda = Combination of (l0/s)*(alpha*Psi)/2
        % Input
        % eta0 = reference viscosity at referenc stress
        % s0   = reference stress
        % l0   = initial lenght
        % D0   = initial thickness
        % n    = power law exponent
        % Df   = viscosity contrast between diffusion and dislocation at reference
        %        stress
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        drho   = 2*s0/(9.81*l0);     % delta rho 
        B_n   = s0^(1-n)/eta0;       % compliance dislocation
        B_d   = 1/(Df*eta0);         % compliance diffusion
        ec    = (B_n*s0^n+B_d*s0);   % characteristic strain rate
        tc    = 1/ec             ;   % characteristic time scale
        etaS_eff = (1/eta0+1/(Df*eta0))^(-1); % effective viscosity of the slab at reference condition
        Psi      = etaum/etaS_eff;    % ratio between the upper mantle viscosity and the slab viscosity
        % Alpha 
        alpha = 5.0;                 % Ancient parameter derived by Yanick et al. 1986
        Len = l0/(2*1000e3);         % Length divided by a characteristic lenght scale (i.e. size of my model)
        Lambda = Len*alpha*Psi;      % Parameter derived by 2D numerical simulation 

end


function plot_results(Tests,name,ptsave)
    % Input 
    % Testdata=> Data Structure containing all the tests 
    % Short description
    % Plot t/tc - Dnorm of all the tests and colored as a function of
    % lambda.
    % Plot td against Lambda scatter plot.
    % Stress => To Do How to retrieve the stress data per each timestep 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % function plot not dimensional data {+ additional collect data for the scatter plot
    % i.e. Lambda, t_d,Psi} 
    [Data_S] =  Plot_1D_Plots(Tests,name,ptsave);
    % function to do scatter plot
    plot_scatter(Data_S,name,ptsave,'t_det');
    plot_scatter(Data_S,name,ptsave,'tau_ii');
    plot_scatter(Data_S,name,ptsave,'time_max-time det');



end

function [Data_S] = Plot_1D_Plots(Tests,name,ptsave)

    % Collect the field names 
    fn = fieldnames(Tests);
    % number of test 
    itest = length(fn);
    % Prepare Data_S array
    Data_S = zeros(6,itest); 
    i = 1;
    for k = 1:numel(fn) 
       
       TD = Tests.(fn{k});
       Data_S(1,i) = TD.initial_data.Lambda;
       Data_S(2,i) = log10(TD.initial_data.Psi); 
       Data_S(3,i) = TD.initial_data.n*TD.t_det;
       Data_S(4,i) = TD.t_t_max;
       Data_S(5,i) = TD.time_t_M*TD.initial_data.n;
       Data_S(6,i) = TD.t_t_det;
       i = i+1;
    end

    tic
    plot1D_setExp(Tests,Data_S,name,'D_norm',ptsave)
    toc
    tic
    plot1D_setExp(Tests,Data_S,name,'tau_eff',ptsave)
    toc
    tic
    plot1D_setExp(Tests,Data_S,name,'tau_D_tau_B',ptsave)
    toc
end

function plot_scatter(Data_S,name,ptsave,field) 
c = Data_S(2,:);
double = 0; 
if strcmp(field,'t_det')
    x = Data_S(1,:);
    y = Data_S(3,:);
    ylabel('t^O_d/t^P_d [n.d]')
    fin = 'Global_test1DS_time_det';
elseif strcmp(field,'tau_ii')
    x = Data_S(1,:);
    y = Data_S(4,:); % tau max
    y2 = Data_S(6,:); % tau @ detachment
    ylabel('$\frac{\tau_{eff}}{\tau_{B_0}}$','Interpreter','latex');
    double = 1.0 ; 
    fin = 'Global_test1DS_Stress';

elseif strcmp(field,'time_max-time det')
    x = Data_S(1,:);
    y = Data_S(3,:)-Data_S(5,:);
    ylabel('$t_{det}-t(\tau_{eff}^{MAX}) [n.d.]$','Interpreter','latex');
    fin = 'Global_test1DS_dT';
end
figure(1)
scatter(x,y,20,c,'filled','d')
hold on 

if double >0 
    scatter(x,y2,5,c,'filled','o')
end

try
    cmap = colormap(crameri('nuuk'));
catch
    cmap = colormap('jet');
end
% Set colorbar
c=colorbar;
c.Label.String = 'log10(\Psi) [n.d.]';
grid on
box on
xlabel('\Lambda [n.d.]')
%ylabel('t^O_d/t^P_d [n.d]')
%ylim([0.8,20])
if strcmp(field,'t_det')
    ylabel('t^O_d/t^P_d [n.d]')
 elseif strcmp(field,'tau_ii')
    ylabel('$\frac{\tau_{eff}}{\tau_{B_0}}$','Interpreter','latex'); 
elseif strcmp(field,'time_max-time det')
    ylabel('$t_{det}-t(\tau_{eff}^{MAX}) [n.d.]$','Interpreter','latex');
 end
set(gca, 'XScale', 'log')
name_picture = strcat(fin,name,'.png');
pt=fullfile(ptsave,'SCAT');
if not(isfolder(pt))
     mkdir(pt);
end
title(name)
pt=fullfile(pt,name_picture);

print(pt,'-dpng')
clf; 
close;

end 
function  plot1D_setExp(Tests, Data_S,name,field,ptsave)

    % Collect the field names 
    fn = fieldnames(Tests);
    % Set the min and max of lambda value for the coloring of the plot
    z_min = log10(min(Data_S(1,:)));
    z_max  =log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
    if strcmp(field,'D_norm')
        VAL = 0;
    elseif strcmp(field,'tau_eff')
        VAL=1; 
    else 
        VAL = 2; 
    end

    try
        cmap = colormap(crameri('Bilbao'));
    catch
        cmap = colormap('jet');
    end
    % Set colorbar
    % Set colorbar
c=colorbar;
c.Label.String = 'log10(\Lambda) [n.d.]';
    % Define coloraxis
    caxis([z_min z_max])
    % min/max are rounded using the log values of Lambda
    Ticks=round((z_min)):round((z_max));
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('10^{%d}',x),Ticks,'UniformOutput',false);
    try
        set(c,'Ticks',Ticks);
        set(c,'TickLabels',TickLabels);
    catch %HG1
        TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
        set(c,'YTick',Ticks);
        set(c,'YTickLabel',TickLabels);
    end
    i = 1; % iterator 
    for k = 1:numel(fn) 
       TD = Tests.(fn{k});
       if VAL == 0 
           buf = TD.D_norm; 
           ylim([0.1,1.0])
           ylabel('$\frac{D}{D_0} [n.d.]$','interpreter','latex')
       elseif VAL == 1
           buf = TD.tau(3,:);
           ylabel('$\frac{\tau_{eff}}{\tau_{B,0}} [n.d.]$','interpreter','latex')
       else
           buf = -TD.tau(2,:)./TD.tau(1,:);
           ylabel('$\frac{\tau_{D}}{\tau_{B}} [n.d.]$','interpreter','latex')
       end

       % Normalize Lambda value w.r.t. the limit that I assumed to be
       % likely
       V = (log10(TD.initial_data.Lambda)-z_min)/(z_max-z_min);
       if V<0 
           V=0;
       elseif V>1
            V=1;
       end
       V=round(1+V*(size(cmap,1)-1));%round to nearest index
       C = cmap(V,:);
       hold on 
       plot(TD.time*TD.initial_data.n,buf,'Color',C)
       grid on 
       xlim([0,20])
       xlabel('$n*(\frac{t}{t_c}) [n.d.]$','interpreter','latex')
       i = i+1; 
    end
    box on
    name_picture = strcat('Global_test1D',field,name,'.png');
    pt=fullfile(ptsave,field);
    if not(isfolder(pt))
         mkdir(pt);
    end
    title(name)

    pt=fullfile(pt,name_picture);

    print(pt,'-dpng')
    clf; 
    close;
end
