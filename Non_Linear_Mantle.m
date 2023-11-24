% Figure Manuscript Non Linear portion
clear all
close all
addpath('Utilities/')
addpath('Utilities/Plot_Class/')
addpath('Utilities/ScientificColourMaps8/')
clf; 
%%NonLinear_Tests_Data_Base_n3_iteration_high_T.mat It=load('NonLinear_Tests_Data_Base_n3_iteration.mat');

It=load('..\Data_Base\NonLinear_Tests_Data_Base_n3_iteration_high_T_initial_guess.mat');
ptsave = '../NonLinearMantle_PP';

%% Are the initial viscosities retrivable from the numerical data? 
S = It.Data_S; 
theta = 80*5.0./1000; 

Lambda0 = S.Lambda0; % Lambda parameter that is derived from the first step numerical step
Lambda  = S.Lambda; 
eta0um  = S.eta_um_0.*S.tc.*S.tau0; 
eps0um  = S.eps_um_0;
tau0um  = S.tau_um_0; 
gamma   = (S.L0.*5.0)./2000e3;
xdisl0  = S.xdisl; 
xium    = S.xiUM; 
n       = S.n; 
Lambda_c = S.Lambda0.*(1+xium.*tau0um.^(n-1));
allowed = ~isnan(S.tdet)==1 & S.xdisl>=0.5;
T       = (1-n)./n;
eta_disl_um = (1./xium).*S.eta0DUM.*(xdisl0.*eps0um).^T;
correction = (1./S.eta0DUM+1./eta_disl_um).^(-1);
correction = correction./S.eta0DUM;
correction2 = 1./((1)+xium.*(xdisl0.*theta).^(-T));
Lambda_c2 = S.Lambda.*correction2;
td = 1./S.tdet; 

Bd_um = S.BdUM.*S.tau0.*S.tc;
Bn_um = S.BnUM.*S.tau0.^S.n.*S.tc;
n     = 3.5; 
[tau_0_m_c,eta_eff_m_c,xdisl_0] = Post_Process_stress_data(Bn_um,xium,S.eta0DUM./S.tau0./S.tc,theta,n);
%%
allowed = ~isnan(S.tdet)==1; %& S.xdisl>=0.5;

figure_lambda0_td = scatter_plot_post_process;
figure_lambda0_td.figure_number = 1; 
figure_lambda0_td.x = S.Lambda0(allowed==1);
figure_lambda0_td.y = 1./S.tdet(allowed==1);
figure_lambda0_td.c = S.xiUM(allowed==1);
figure_lambda0_td.colormap_discrete = 8; 
figure_lambda0_td.logx = 'log';
figure_lambda0_td.logy = 'linear';
figure_lambda0_td.logcolor = 1; 
figure_lambda0_td.colormap_f = 'oslo';
figure_lambda0_td.xlabel = '$\Lambda_0$';
figure_lambda0_td.ylabel = '$t^{\dagger,-1}_d$';
figure_lambda0_td.clabel = '${\nu_{M,0}}$';
figure_lambda0_td.save_path = ptsave;% ptsave;
figure_lambda0_td.name_figure = 'Lambda_0_vs_td';
figure_lambda0_td.size_picture = [12.5,13];
figure_lambda0_td.make_scatter_plot; 
%%
allowed = ~isnan(S.tdet)==1; %& S.xdisl>=0.5;

figure_lambdac_td = scatter_plot_post_process;
figure_lambdac_td.figure_number = 1; 
figure_lambdac_td.x = S.Lambda(allowed==1);
figure_lambdac_td.y = 1./S.tdet(allowed==1);
figure_lambdac_td.c = S.xiUM(allowed==1);
figure_lambdac_td.colormap_discrete = 10; 
figure_lambdac_td.logx = 'log';
figure_lambdac_td.logy = 'linear';
figure_lambdac_td.logcolor = 0; 
figure_lambdac_td.colormap_f = 'glasgow';
figure_lambdac_td.xlabel = '$\Lambda_c$';
figure_lambdac_td.ylabel = '$t^{\dagger,-1}_d$';
figure_lambdac_td.clabel = '${\nu_{M,0}}$';
figure_lambdac_td.save_path = ptsave;% ptsave;
figure_lambdac_td.name_figure = 'Lambda_0_vs_td';
figure_lambdac_td.size_picture = [12.5,13];
figure_lambdac_td.make_scatter_plot; 

%%
allowed = ~isnan(S.tdet)==1; %& S.xdisl>=0.5;
Psi = (1+S.xiUS).*S.eta0DUM./(S.eta0DS);
cor = 1+S.xiUM.*S.tau_um_0.^(S.n-1);
Psi_cor = Psi./cor; 

Psi_C_fig = scatter_plot_post_process;
Psi_C_fig.figure_number = 1; 
Psi_C_fig.x = Psi(allowed==1);
Psi_C_fig.y = S.tau_um_0(allowed==1);
Psi_C_fig.c = S.Lambda0(allowed==1);
Psi_C_fig.logx = 'log';
Psi_C_fig.logy = ['log'];
Psi_C_fig.logcolor = 1;
Psi_C_fig.xlim=[10^(-9),5*10^5];
cbar1 = log10(S.xiUM);
%Psi_C_fig.colormap_discrete = length(min(cbar1):1:max(cbar1))-1;
Psi_C_fig.clim = [-4,0]
Psi_C_fig.colormap_f = 'lipari';
Psi_C_fig.xlabel = '$\Psi_c$';
Psi_C_fig.ylabel = '$\tau_{M.0}$';
Psi_C_fig.clabel = '${\xi_{M}}$';
Psi_C_fig.save_path = ptsave;% '../Post_process2/';
Psi_C_fig.name_figure = 'Psi_cvstauM';
Psi_C_fig.size_picture = [12.5,13];
Psi_C_fig.make_scatter_plot; 

%%
allowed = ~isnan(S.tdet)==1; %& S.xdisl>=0.5;
Psi = (1+S.xiUS).*S.eta0DUM./(S.eta0DS);
cor = 1+S.xiUM.*S.tau_um_0.^(S.n-1);
Psi_cor = Psi./cor; 

Psi_0_fig = scatter_plot_post_process;
Psi_0_fig.figure_number = 1; 
Psi_0_fig.x = Psi_cor(allowed==1);
Psi_0_fig.y = S.tau_um_0(allowed==1);
Psi_0_fig.c = S.Lambda0(allowed==1);
Psi_0_fig.logx = 'log';
Psi_0_fig.logy = ['log'];
Psi_0_fig.logcolor = 1;
cbar1 = log10(S.xiUM);
Psi_0_fig.colormap_discrete = length(min(cbar1):1:max(cbar1))-1;
Psi_0_fig.xlim=[10^(-9),5*10];

Psi_0_fig.colormap_f = 'lipari';
Psi_0_fig.xlabel = '$\Psi_c$';
Psi_0_fig.ylabel = '$\tau_{M.0}$';
Psi_0_fig.clabel = '${\xi_{M}}$';
Psi_0_fig.save_path = ptsave;% '../Post_process2/';
Psi_0_fig.name_figure = 'Psi_0vstauM';
Psi_0_fig.size_picture = [12.5,13];
Psi_0_fig.make_scatter_plot; 


%%
allowed = ~isnan(S.tdet)==1;% & (S.xdisl)<0.5;
Psi = (1+S.xiUS).*S.eta0DUM./(S.eta0DS);
cor = 1+S.xiUM.*S.tau_um_0.^(S.n-1);
Psi_cor = Psi./cor; 

figure_lambda0_lambda_3 = scatter_plot_post_process;
figure_lambda0_lambda_3.figure_number = 1; 
figure_lambda0_lambda_3.x = (S.Lambda0(allowed==1));
figure_lambda0_lambda_3.y = (S.tau_drag_initial(allowed==1));
figure_lambda0_lambda_3.c = (S.xdisl(allowed==1));
figure_lambda0_lambda_3.logx = 'log';
figure_lambda0_lambda_3.logy = ['linear'];
figure_lambda0_lambda_3.logcolor = 0; 
figure_lambda0_lambda_3.colormap_f = 'glasgow';
figure_lambda0_lambda_3.xlabel = '$\Lambda_0$';
figure_lambda0_lambda_3.ylabel = '$\tau_{D,0}$';
figure_lambda0_lambda_3.clabel = '${\xi_{M}}$';
figure_lambda0_lambda_3.save_path = ptsave;% '../Post_process/';
figure_lambda0_lambda_3.name_figure = 'noname';
figure_lambda0_lambda_3.size_picture = [12.5,13];
figure_lambda0_lambda_3.make_scatter_plot; 

%%
Lambda_3 = S.Lambda./(1+xium.*tau0um.^(S.n-1));
figure_lambda0_lambda_3 = scatter_plot_post_process;
figure_lambda0_lambda_3.figure_number = 12; 
figure_lambda0_lambda_3.x = S.Lambda0(allowed==1);
figure_lambda0_lambda_3.y = 1./S.tdet(allowed==1);
figure_lambda0_lambda_3.c = log10(S.xiUM(allowed==1));
figure_lambda0_lambda_3.logx = 'log';
figure_lambda0_lambda_3.logy = 'linear';
figure_lambda0_lambda_3.logcolor = 0; 
figure_lambda0_lambda_3.colormap_f = 'glasgow';
figure_lambda0_lambda_3.xlabel = '$\frac{\Lambda_c}{1+\xi_M\tau^{\dagger,n-1}_{M,0}}$';
figure_lambda0_lambda_3.ylabel = '$\Lambda_0$';
figure_lambda0_lambda_3.clabel = '${\xi_{M}}$';
figure_lambda0_lambda_3.save_path = ptsave;% '../Post_process/';
figure_lambda0_lambda_3.name_figure = 'Lambda_c_vs_Lambda_0_c_corrected';
figure_lambda0_lambda_3.size_picture = [12.5,13];
figure_lambda0_lambda_3.make_scatter_plot; 
%%
allowed = ~isnan(S.tdet)==1 ;%& (S.xiUM)==1e-2;
Lambda_3 = S.Lambda./(1+xium.*tau0um.^(S.n-1));
figure_lambda0_lambda_3 = scatter_plot_post_process;
figure_lambda0_lambda_3.figure_number = 1; 
figure_lambda0_lambda_3.x = Lambda0(allowed==1);
figure_lambda0_lambda_3.y = 1./S.tdet(allowed==1);
figure_lambda0_lambda_3.c = S.xiUM(allowed==1);
cbar1 = log10(S.xiUM);
figure_lambda0_lambda_3.colormap_discrete = length(min(cbar1):1:max(cbar1))-1;
figure_lambda0_lambda_3.xlim = [10^(-8),5.0];
figure_lambda0_lambda_3.ylim = [0.0,1.1];
figure_lambda0_lambda_3.logx = 'log';
figure_lambda0_lambda_3.logy = 'linear';
figure_lambda0_lambda_3.logcolor = 1; 
figure_lambda0_lambda_3.colormap_f = 'managua';
figure_lambda0_lambda_3.xlabel = '$\Lambda_0$';
figure_lambda0_lambda_3.ylabel = '$\frac{1}{t_d}$';
figure_lambda0_lambda_3.clabel = '$log_{10}{\Lambda_c}$';
figure_lambda0_lambda_3.save_path = ptsave;% '../Post_process/';
figure_lambda0_lambda_3.name_figure = 'Lambda_0_vs_td';
figure_lambda0_lambda_3.size_picture = [12.5,13];
figure_lambda0_lambda_3.make_scatter_plot; 
%%
figure_lambda_3_td = scatter_plot_post_process;
figure_lambda_3_td.figure_number = 2; 
figure_lambda_3_td.x = S.Lambda0(allowed==1);
figure_lambda_3_td.y = 1./td(allowed==1);
figure_lambda_3_td.c = xdisl_0(allowed==1);
figure_lambda_3_td.logx = 'log';
figure_lambda_3_td.logy = 'linear';
figure_lambda_3_td.logcolor = 0; 
figure_lambda_3_td.colormap_f = 'glasgow';
figure_lambda_3_td.xlabel = '$\Lambda_c$';
figure_lambda_3_td.ylabel = '$\frac{1}{t^{\dagger}_d} $';
figure_lambda_3_td.clabel = '$log_{10}{\nu_c^M}$';
figure_lambda_3_td.save_path = ptsave;% '../Post_process/';
figure_lambda_3_td.name_figure = 'detachment_timescale';
figure_lambda_3_td.size_picture = [12.5,13];
figure_lambda_3_td.make_scatter_plot; 
%%
figure_lambda_3_tau_max = scatter_plot_post_process;
figure_lambda_3_tau_max.figure_number = 3; 
figure_lambda_3_tau_max.x = S.Lambda0(allowed==1);
figure_lambda_3_tau_max.y = S.tau_max(allowed==1);
figure_lambda_3_tau_max.c = S.xdisl(allowed==1);
figure_lambda_3_tau_max.xlim = [10^(-8),5.0];
figure_lambda_3_tau_max.ylim = [0.0,10.5];
figure_lambda_3_tau_max.logx = 'log';
figure_lambda_3_tau_max.logy = 'linear';
figure_lambda_3_tau_max.logcolor = 0; 
figure_lambda_3_tau_max.colormap_f = 'glasgow';
figure_lambda_3_tau_max.xlabel = '$\Lambda_c$';
figure_lambda_3_tau_max.ylabel = '$\tau^{Max,\dagger}$';
figure_lambda_3_tau_max.clabel = '$log_{10}{\xi^{M}}$';
figure_lambda_3_tau_max.save_path = ptsave;% '../Post_process/';
figure_lambda_3_tau_max.name_figure = 'tau_max';
figure_lambda_3_tau_max.size_picture = [12.5,13];
figure_lambda_3_tau_max.make_scatter_plot; 

%% 
eta0D = 10.^(-7:0.1:3);
xium  = 10.^(-10:0.1:10);
[eta,Xi] = meshgrid(eta0D,xium); 
theta = 80*5.0./1000; 
xdisl_art = 0.0.*eta; 
tau_m     = 0.0.*eta; 
eta_mantle = 0.0.*eta; 

n     = 3.5; 
Gamma  = (1-n)./n;
options = optimset('TolX',1e-6);

for i = 1:length(Xi(:))
    A = cputime;
    % Find the stress using fzero    
    % Since I am effectively mentally retarded, i was using my usual
    % approach for finding the dislocation creep rheology, neglecting the
    % small and pesky details that the characteristic strain rate is
    % defined using a reference stress, so, the strain rate of the upper
    % mantle is not the same of the one that you might find in the slab. 
    % Computing the bn 
    Bn = 1./(2.*eta(i)./Xi(i));
    % Computing bd
    Bd = 1./(2.*eta(i));
    % Computing the strain rate indipendent part
    eta0=0.5 * Bn^(-1/n);
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
    eps_disl = Bn.*tau_m(i).^n;
    %Computing the partition
    xdisl_art(i) = eps_disl./theta;
    %Computing the real effective viscosity
    eta_mantle(i) = 2.*theta.*tau_m(i);
    B = cputime;
    disp(['time =',num2str(B-A), 's']);
end

%%
Figure_partition = maps_plot;
Figure_partition.figure_number = 4; 
Figure_partition.contourf_option={1,10};
Figure_partition.x = eta;
Figure_partition.y = Xi; 
Figure_partition.c = xdisl_art; 
Figure_partition.colormap_f = 'managua';
Figure_partition.logx = 'log';
Figure_partition.logy = 'log';
Figure_partition.logcolor = 0; 
Figure_partition.xlabel = '$\eta^{\dagger}_{D}$';
Figure_partition.ylabel = '$\xi$';
Figure_partition.clabel = '$\nu_{disl}$';
Figure_partition.name_figure = 'partition_map_plot';
Figure_partition.save_path = ptsave;% '../Post_process/';
Figure_partition.clim      = [0,1.0];
Figure_partition.size_picture = [12.5,13];
Figure_partition.make_maps;

Figure_viscosity = maps_plot;
Figure_viscosity.figure_number = 5; 
Figure_viscosity.contourf_option={1,10};
Figure_viscosity.x = eta;
Figure_viscosity.y = Xi; 
Figure_viscosity.c = eta_mantle./eta; 
Figure_viscosity.colormap_f = 'navia';
Figure_viscosity.logx = 'log';
Figure_viscosity.logy = 'log';
Figure_viscosity.logcolor = 1; 
Figure_viscosity.xlabel = '$\eta^{\dagger}_{D}$';
Figure_viscosity.ylabel = '$\xi$';
Figure_viscosity.clabel = '$log_{10}\left(\frac{\eta^{\dagger}_{eff}}{\eta^{\dagger}_{D}}\right)$';
Figure_viscosity.name_figure = 'viscosity_map_plot';
Figure_viscosity.save_path = ptsave;% '../Post_process/';
lim=prctile(eta_mantle,[10 90],'all');
Figure_viscosity.clim      = [log10(lim')];
Figure_viscosity.size_picture = [12.5,13];
Figure_viscosity.make_maps;



Figure_stress = maps_plot;
Figure_stress.figure_number = 6; 
Figure_stress.contourf_option={1,10};
Figure_stress.x = eta;
Figure_stress.y = Xi; 
Figure_stress.c = tau_m; 
Figure_stress.colormap_f = 'lipari';
Figure_stress.logx = 'log';
Figure_stress.logy = 'log';
Figure_stress.logcolor = 1; 
Figure_stress.xlabel = '$\eta^{\dagger}_{D}$';
Figure_stress.ylabel = '$\xi$';
Figure_stress.clabel = '$log_{10}\left(\tau\right)$';
Figure_stress.name_figure = 'tau_map';
Figure_stress.save_path = ptsave;% '../Post_process/';
lim=prctile(tau_m,[10 90],'all');
Figure_stress.clim      = [log10(lim')];
Figure_stress.size_picture = [12.5,13];
Figure_stress.make_maps;

%%
C = (eta./(1+Xi.*tau_m.^(n-1))).*0.75;

Figure_stress = maps_plot;
Figure_stress.figure_number = 6; 
Figure_stress.contourf_option={0,10};
Figure_stress.x = eta;
Figure_stress.y = Xi; 
Figure_stress.c = C; 
Figure_stress.colormap_f = 'lipari';
Figure_stress.logx = 'log';
Figure_stress.logy = 'log';
Figure_stress.logcolor = 1; 
Figure_stress.xlabel = '$\eta^{\dagger}_{D}$';
Figure_stress.ylabel = '$\xi$';
Figure_stress.clabel = '$log_{10}\left(\frac{\eta_0 \gamma}{1+\xi\tau_{M,c}.^{n-1}}\right)$';
Figure_stress.name_figure = 'Lambda_map';
Figure_stress.save_path = ptsave;% '../Post_process/';
lim=prctile(tau_m,[10 90],'all');
Figure_stress.clim      = [log10(lim')];
Figure_stress.size_picture = [12.5,13];
Figure_stress.make_maps;

%% Line Plot 
fun_0D = Manuscript_function_Container;  
% Function that creates the data structure with the relevant information
% for the picture at hand 
path2colormap = strcat('Utilities\ScientificColourMaps8\','lipari','\','lipari','.mat');

load(path2colormap);

cmap        = colormap(lipari);

[Fig1_A]=fun_0D.select_tests_prepare_variables(It.n_3,0,'time_nd','etaum','Lambda0','NonLinear',3,'xius',tau_0_m_c);

[Fig1_B]=fun_0D.select_tests_prepare_variables(It.n_3,0,'time_nd','tau_eff','Lambda0','NonLinear',3,'xius',tau_0_m_c);

[Fig1_C]=fun_0D.select_tests_prepare_variables(It.n_3,0,'time_nd','D_norm','Lambda0','NonLinear',3,'xius',tau_0_m_c);

[Fig1_D]=fun_0D.select_tests_prepare_variables(It.n_3,0,'time_nd','xdisl','xium','NonLinear',3,'xius',tau_0_m_c);

[Fig1_E]=fun_0D.select_tests_prepare_variables(It.n_3,0,'time_nd','tau_M','xium','NonLinear',3,'xius',tau_0_m_c);
%%

size_tests = length(squeeze(Fig1_E(1,1,:)));

c = squeeze(Fig1_E(3,:,:));

color_lists = fun_0D.color_computation(size_tests,c,-4,0);
Figure_viscosity_tests = line_plot_post_process;
Figure_viscosity_tests.figure_number = 7; 
Figure_viscosity_tests.logx = 'linear';
Figure_viscosity_tests.logy = 'log';
Figure_viscosity_tests.colormap_f = cmap; 
Figure_viscosity_tests.xlabel = '$t^{\dagger}$';
Figure_viscosity_tests.ylabel = '$\tau_{M}$';
Figure_viscosity_tests.clabel = '$log_{10}\left(\Lambda_c\right)$';
Figure_viscosity_tests.clim  = [-4,0];
Figure_viscosity_tests.ctick = [-4:1:0];
Figure_viscosity_tests.xlim  = [0.01,40];
Figure_viscosity_tests.x      = squeeze(Fig1_C(1,:,:));
Figure_viscosity_tests.y      = squeeze(Fig1_C(2,:,:));
Figure_viscosity_tests.c      = log10(squeeze(Fig1_C(3,:,:)));
Figure_viscosity_tests.size_picture = [13,12.5];
Figure_viscosity_tests.name_figure = 'tau_Mevolution';
Figure_viscosity_tests.save_path = ptsave;% '../Post_process/';

Figure_viscosity_tests.multiline_plot(size_tests,color_lists)

%%

size_tests = length(squeeze(Fig1_B(1,1,:)));

c = squeeze(Fig1_B(3,:,:));

color_lists = fun_0D.color_computation(size_tests,c,-4,0);
Figure_viscosity_tests = line_plot_post_process;
Figure_viscosity_tests.figure_number = 7; 
Figure_viscosity_tests.logx = 'log';
Figure_viscosity_tests.logy = 'log';
Figure_viscosity_tests.colormap_f = cmap; 
Figure_viscosity_tests.xlabel = '$t^{\dagger}$';
Figure_viscosity_tests.ylabel = '${\eta_{M}^{\dagger}}$';
Figure_viscosity_tests.clabel = '$log_{10}\left(\Lambda_c\right)$';
Figure_viscosity_tests.clim  = [-4,0];
Figure_viscosity_tests.ctick = [-4:1:0];
Figure_viscosity_tests.xlim  = [0.01,40];
Figure_viscosity_tests.ylim  = [10^-1,30];
Figure_viscosity_tests.x      = squeeze(Fig1_A(1,:,:));
Figure_viscosity_tests.y      = squeeze(Fig1_A(2,:,:));
Figure_viscosity_tests.c      = log10(squeeze(Fig1_A(3,:,:)));
Figure_viscosity_tests.size_picture = [13,12.5];
Figure_viscosity_tests.name_figure = 'Viscosity_evolution';
Figure_viscosity_tests.save_path = ptsave;% '../Post_process/';

Figure_viscosity_tests.multiline_plot(size_tests,color_lists)
%%
Figure_stress_tests = line_plot_post_process;
Figure_stress_tests.figure_number = 8; 
Figure_stress_tests.logx = 'log';
Figure_stress_tests.logy = 'log';
Figure_stress_tests.colormap_f = cmap; 
Figure_stress_tests.xlabel = '$t^{\dagger}$';
Figure_stress_tests.ylabel = '$\tau^{\dagger}$';
Figure_stress_tests.clabel = '$log_{10}\left(\Lambda_0\right)$';
Figure_stress_tests.clim  = [-4,0];
Figure_stress_tests.ctick = [-4:1:0];
Figure_stress_tests.xlim  = [0.01,40];
Figure_stress_tests.x      = squeeze(Fig1_B(1,:,:));
Figure_stress_tests.y      = squeeze(Fig1_B(2,:,:));
Figure_stress_tests.c      = log10(squeeze(Fig1_B(3,:,:)));
Figure_stress_tests.size_picture = [13,12.5];
Figure_stress_tests.name_figure = 'Stress_evolution';
Figure_stress_tests.save_path = ptsave;% '../Post_process/';

Figure_stress_tests.multiline_plot(size_tests,color_lists)

%%

nu = squeeze(Fig1_C(4,1,:));
x  = squeeze(Fig1_C(1,:,:));
y  = squeeze(Fig1_C(2,:,:));
z = squeeze(Fig1_C(3,:,:));
xdf = x(:,nu<0.4);
ydf = y(:,nu<0.4);
zdf = z(:,nu<0.4);

color_lists = fun_0D.color_computation(length(nu(nu<0.4)),zdf,-4,0);


Figure_thickness = line_plot_post_process;
Figure_thickness.figure_number = 9; 
Figure_thickness.logx = 'linear';
Figure_thickness.logy = 'linear';
Figure_thickness.colormap_f = cmap; 
Figure_thickness.xlabel = '$t^{\dagger}$';
Figure_thickness.ylabel = '$D^{\dagger}$';
Figure_thickness.clabel = '$log_{10}\left(\Lambda_0\right)$';
Figure_thickness.clim  = [-4,0];
Figure_thickness.ctick = [-4:1:0];
Figure_thickness.xlim  = [0.0,40];
Figure_thickness.ylim  = [0.1,1.1];
Figure_thickness.x      = xdf;
Figure_thickness.y      = ydf;
Figure_thickness.c      = log10(zdf);
Figure_thickness.size_picture = [13,12.5];
Figure_thickness.name_figure = 'Thickness_evolution_DIFFUSION';
Figure_thickness.save_path = ptsave;% '../Post_process/';
Figure_thickness.multiline_plot(length(nu(nu<0.4)),color_lists)

%%
nu = squeeze(Fig1_E(4,1,:));
x  = squeeze(Fig1_E(1,:,:));
y  = squeeze(Fig1_E(2,:,:));
z = squeeze(Fig1_E(3,:,:));
xdf = x(:,nu>=0.5);
ydf = y(:,nu>=0.5);
zdf = z(:,nu>=0.5);

color_lists = fun_0D.color_computation(length(nu(nu>=0.4)),zdf,-2,6);

path2colormap = strcat('Utilities\ScientificColourMaps8\','batlowK','\','batlowK','.mat');

load(path2colormap);

cmap        = colormap(batlowK);
Fstress = line_plot_post_process;
Fstress.figure_number = 10; 
Fstress.logx = 'linear';
Fstress.logy = 'log';
Fstress.colormap_f = cmap; 
Fstress.xlabel = '$t^{\dagger}$';
Fstress.ylabel = '$\tau_{M}^{\dagger}$';
Fstress.clabel = '$log_{10}\left(\xi_M\right)$';
Fstress.clim  = [-2,6];
Fstress.ctick = [-2:1:6];
Fstress.xlim  = [0.0,20];
%Fstress.ylim  = [0.1,1.1];
Fstress.x      = xdf;
Fstress.y      = ydf;
Fstress.c      = log10(zdf);
Fstress.size_picture = [13,12.5];
Fstress.name_figure = 'stress_evolution_dislocation';
Fstress.save_path = ptsave;% '../Post_process/';
Fstress.multiline_plot(length(nu(nu>=0.4)),color_lists)

%%
nu = squeeze(Fig1_E(4,1,:));
x  = squeeze(Fig1_E(1,:,:));
y  = squeeze(Fig1_E(2,:,:));
z = squeeze(Fig1_E(3,:,:));
xdf = x(:,nu<0.5);
ydf = y(:,nu<0.5);
zdf = z(:,nu<0.5);

color_lists = fun_0D.color_computation(length(nu(nu>=0.4)),zdf,-2,6);

path2colormap = strcat('Utilities\ScientificColourMaps8\','batlowK','\','batlowK','.mat');

load(path2colormap);

cmap        = colormap(batlowK);
Fstress = line_plot_post_process;
Fstress.figure_number = 10; 
Fstress.logx = 'linear';
Fstress.logy = 'log';
Fstress.colormap_f = cmap; 
Fstress.xlabel = '$t^{\dagger}$';
Fstress.ylabel = '$\tau_{M}^{\dagger}$';
Fstress.clabel = '$log_{10}\left(\xi_M\right)$';
Fstress.clim  = [-2,6];
Fstress.ctick = [-2:1:6];
Fstress.xlim  = [0.0,20];
%Fstress.ylim  = [0.1,1.1];
Fstress.x      = xdf;
Fstress.y      = ydf;
Fstress.c      = log10(zdf);
Fstress.size_picture = [13,12.5];
Fstress.name_figure = 'stress_evolution_diffusion';
Fstress.save_path = ptsave;% '../Post_process/';
Fstress.multiline_plot(length(nu(nu>=0.4)),color_lists)

%%

%%
x = squeeze(Fig1_D(1,:,:));
y = squeeze(Fig1_D(2,:,:));
z = squeeze(Fig1_D(3,:,:));
c = squeeze(Fig1_D(4,:,:));
tt = squeeze(Fig1_B(2,:,:));
Lambda = squeeze(Fig1_C(3,:,:));
lambda = zeros(size_tests,1);
equal_td = zeros(size_tests,1);
max_     = zeros(size_tests,1);
L         = zeros(size_tests,1);
td        = zeros(size_tests,1); 
tau_max   = zeros(size_tests,1);
cc = zeros(size_tests,1);
for i=1:size_tests
    ind = find(y(:,i)==nanmax(y(:,i)),1);
    if ~isempty(ind)
        equal_td(i) = nanmean(y(:,i));
        max_(i)      = nanmax(y(:,i));
        L(i)        = c(1,i);
        cc(i)       = z(1,i);
        ind_time = find(isnan(y(:,i)),1)-1;
        td(i)       = x(ind_time,i);
        tau_max(i)  =  max(tt(:,i));
        lambda(i)    = Lambda(1,i);
    else
        equal_td = NaN; 
        L(i) = NaN; 
        cc(i) = NaN;
        max_(i) = Nan; 
        td(i)       = NaN;
        tau_max(i)  =  NaN;
        lambda(i) = NaN; 
    end
end



cbar1 = log10(S.xiUM);

fig_s = scatter_plot_post_process;
fig_s.figure_number = 3; 
fig_s.markers = {'d','o'};
fig_s.conditions = {max_,0.5,'leq'}; 
fig_s.x = L;
fig_s.y = equal_td-0.5;
fig_s.c = cc;
fig_s.xlim = [min(L)-0.5*min(L),max(L)+0.5*max(L)];
fig_s.ylim = [-0.55,0.55];
fig_s.colormap_discrete = length(min(cbar1):1:max(cbar1))-1;
fig_s.logx = 'log';
fig_s.logy = 'linear';
fig_s.logcolor = 1; 
fig_s.colormap_f = 'managua';
fig_s.xlabel = '$\Lambda_c$';
fig_s.ylabel = '$\tau^{Max,\dagger}$';
fig_s.clabel = '$log_{10}{\xi^{M}}$';
fig_s.save_path = ptsave;% '../Post_process/';
fig_s.name_figure = 'tau_max';
fig_s.clim  = [-2,6];
fig_s.size_picture = [12.5,13];
fig_s.make_scatter_plot; 
figure(3)
ax = gca; 
yline(0.0,'r','LineWidth',2.0)
xline(1.0,'g','LineWidth',2.0,'Alpha',0.2,'LineStyle',':')
filename = fullfile(fig_s.save_path);

filename = fullfile(filename,strcat(fig_s.name_figure,'.png'));
exportgraphics(ax,filename,'Resolution',600,'BackgroundColor','white') 


path2colormap = strcat('Utilities\ScientificColourMaps8\','glasgow','\','glasgow','.mat');

load(path2colormap);
%%
cmap        = colormap(lipari);

size_tests = length(squeeze(Fig1_E(1,1,:)));

c = squeeze(Fig1_D(3,:,:));

allowed = max_<=0.5 & equal_td <= 0.5; 
x = squeeze(Fig1_D(1,:,:));
y =  squeeze(Fig1_D(2,:,:));
c = squeeze(Fig1_C(3,:,:));

color_lists = fun_0D.color_computation(length(allowed(allowed==1)),c(:,allowed==1),-4,0);

% Effectively linear experiments

Figure_partition = line_plot_post_process;
Figure_partition.figure_number = 11; 
Figure_partition.logx = 'linear';
Figure_partition.logy = 'linear';
Figure_partition.colormap_f = cmap; 
Figure_partition.xlabel = '$t^{\dagger}$';
Figure_partition.ylabel = '$D^{\dagger}$';
Figure_partition.clabel = '$log_{10}\left(\Lambda_0\right)$';
Figure_partition.clim  = [-4,0];
Figure_partition.ctick = [-4:1:0];
Figure_partition.xlim  = [0.01,40];
Figure_partition.ylim  = [0.0,1.0];
Figure_partition.x      = x(:,allowed==1);
Figure_partition.y      = y(:,allowed==1);
Figure_partition.c      = log10(c(:,allowed==1));
Figure_partition.size_picture = [13,12.5];
Figure_partition.name_figure = 'nu_Linear';
Figure_partition.save_path = ptsave;% '../Post_process/';

Figure_partition.multiline_plot(length(allowed(allowed==1)),color_lists)
figure(11)
ax = gca; 

yline(0.5,'LineWidth',2.0,'alpha',0.5);
filename = fullfile(Figure_partition.save_path);
if ~isdir(filename)
    mkdir(filename);
end
filename = fullfile(filename,strcat(Figure_partition.name_figure,'.png'));
exportgraphics(ax,filename,'Resolution',600,'BackgroundColor','white')


% Mixed Cases

allowed = max_>0.5 & equal_td <= 0.5; 
x = squeeze(Fig1_D(1,:,:));
y =  squeeze(Fig1_D(2,:,:));
c = squeeze(Fig1_C(3,:,:));

color_lists = fun_0D.color_computation(length(allowed(allowed==1)),c(:,allowed==1),-4,0);

% Effectively linear experiments

Figure_partition = line_plot_post_process;
Figure_partition.figure_number = 12; 
Figure_partition.logx = 'linear';
Figure_partition.logy = 'linear';
Figure_partition.colormap_f = cmap; 
Figure_partition.xlabel = '$t^{\dagger}$';
Figure_partition.ylabel = '$D^{\dagger}$';
Figure_partition.clabel = '$log_{10}\left(\Lambda_0\right)$';
Figure_partition.clim  = [-4,0];
Figure_partition.ctick = [-4:1:0];
Figure_partition.xlim  = [0.01,40];
Figure_partition.ylim  = [0.0,1.0];
Figure_partition.x      = x(:,allowed==1);
Figure_partition.y      = y(:,allowed==1);
Figure_partition.c      = log10(c(:,allowed==1));
Figure_partition.size_picture = [13,12.5];
Figure_partition.name_figure = 'nu_Mixed_L';
Figure_partition.save_path = ptsave;% '../Post_process/';

Figure_partition.multiline_plot(length(allowed(allowed==1)),color_lists)
figure(12)
ax = gca; 

yline(0.5,'LineWidth',2.0,'alpha',0.5);
filename = fullfile(Figure_partition.save_path);
if ~isdir(filename)
    mkdir(filename);
end
filename = fullfile(filename,strcat(Figure_partition.name_figure,'.png'));
exportgraphics(ax,filename,'Resolution',600,'BackgroundColor','white')


% Fully non linear

allowed = max_>0.5 & equal_td > 0.5; 
x = squeeze(Fig1_D(1,:,:));
y =  squeeze(Fig1_D(2,:,:));
c = squeeze(Fig1_C(3,:,:));

color_lists = fun_0D.color_computation(length(allowed(allowed==1)),c(:,allowed==1),-4,0);

% Effectively linear experiments

Figure_partition = line_plot_post_process;
Figure_partition.figure_number = 13; 
Figure_partition.logx = 'linear';
Figure_partition.logy = 'linear';
Figure_partition.colormap_f = cmap; 
Figure_partition.xlabel = '$t^{\dagger}$';
Figure_partition.ylabel = '$D^{\dagger}$';
Figure_partition.clabel = '$log_{10}\left(\Lambda_0\right)$';
Figure_partition.clim  = [-4,0];
Figure_partition.ctick = [-4:1:0];
Figure_partition.xlim  = [0.01,40];
Figure_partition.ylim  = [0.0,1.0];
Figure_partition.x      = x(:,allowed==1);
Figure_partition.y      = y(:,allowed==1);
Figure_partition.c      = log10(c(:,allowed==1));
Figure_partition.size_picture = [13,12.5];
Figure_partition.name_figure = 'nu_Non_Linear';
Figure_partition.save_path = ptsave;% '../Post_process/';

Figure_partition.multiline_plot(length(allowed(allowed==1)),color_lists)
figure(13)
ax = gca; 
yline(0.5,'LineWidth',2.0,'alpha',0.5);
filename = fullfile(Figure_partition.save_path);
if ~isdir(filename)
    mkdir(filename);
end
filename = fullfile(filename,strcat(Figure_partition.name_figure,'.png'));
exportgraphics(ax,filename,'Resolution',600,'BackgroundColor','white')


%%
figure(1000)
set(gcf, 'Units','centimeters', 'Position', [0, 0, 13,12.5], 'PaperUnits', 'centimeters', 'PaperSize', [13, 12.5])
clf; 
ax = gca; 
allowed = max_>0.5 & equal_td > 0.5; 
scatter((lambda(allowed==1)),1./td(allowed==1),70,[255,248,220]./255,"filled","square",'MarkerEdgeColor','k')
hold on
allowed = max_>0.5 & equal_td <= 0.5; 
scatter((lambda(allowed==1)),1./td(allowed==1),50,[205,105,201]./255,"filled","d",'MarkerEdgeColor','k')
allowed = max_<=0.5 & equal_td <= 0.5; 
scatter((lambda(allowed==1)),1./td(allowed==1),30,[100,149,237]./255,"filled","o",'MarkerEdgeColor','k')
ax.XScale = 'log'; 
l = legend('Effectively nonlinear','Mixed','Effectively linear','Location','southwest');
l.Interpreter = 'latex';
ax.Box = 'on';
ax.LineWidth = 1.2; 
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.YLim   = [-0.01,1.1];

filename = fullfile(fig_s.save_path);

filename = fullfile(filename,strcat('td','.png'));
exportgraphics(ax,filename,'Resolution',600,'BackgroundColor','white')


figure(101)
set(gcf, 'Units','centimeters', 'Position', [0, 0, 13,12.5], 'PaperUnits', 'centimeters', 'PaperSize', [13, 12.5])
clf; 
ax = gca; 
allowed = max_>0.5 & equal_td > 0.5; 
scatter((lambda(allowed==1)),tau_max(allowed==1),70,[255,248,220]./255,"filled","square",'MarkerEdgeColor','k')
hold on
allowed = max_>0.5 & equal_td <= 0.5; 
scatter((lambda(allowed==1)),tau_max(allowed==1),50,[205,105,201]./255,"filled","d",'MarkerEdgeColor','k')
allowed = max_<=0.5 & equal_td <= 0.5; 
scatter((lambda(allowed==1)),tau_max(allowed==1),30,[100,149,237]./255,"filled","o",'MarkerEdgeColor','k')
ax.XScale = 'log'; 
ax.Box = 'on';
ax.LineWidth = 1.2; 
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
filename = fullfile(fig_s.save_path);

filename = fullfile(filename,strcat('tau_max','.png'));
exportgraphics(ax,filename,'Resolution',600,'BackgroundColor','white') 
%ax.YLim   = [-0.01,1.1];


