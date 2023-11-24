% Figure Realistic case
%% Reshape the datainto matrices
clear all; 
close all; 

% add the path to the main function
addpath('Utilities/')
addpath('Utilities/Plot_Class/')
addpath('Utilities/ScientificColourMaps8/')
folder = '../Manuscript_figure_folder/';
ptsave = '../Realistic_Cases/';
% General information of the size of each figure
if ~isfolder(ptsave)
    mkdir(ptsave);
end
load('../Data_Base/REAL_DATA_tp2');
Vn_test = reshape(Meta_data_Set_Tests.Vnv,Meta_data_Set_Tests.Shape);
Vd_test = reshape(Meta_data_Set_Tests.Vdv,Meta_data_Set_Tests.Shape);
TS_test = reshape(Meta_data_Set_Tests.TSlab,Meta_data_Set_Tests.Shape);
L0_tests      = reshape(Meta_data_Set_Tests.L0,Meta_data_Set_Tests.Shape);
suc_test     = reshape(Meta_data_Set_Tests.suc,Meta_data_Set_Tests.Shape);
td_test     = reshape(Meta_data_Set_Tests.td,Meta_data_Set_Tests.Shape);
Lambda_test     = reshape(Meta_data_Set_Tests.Lambda,Meta_data_Set_Tests.Shape);
tc_test     = reshape(Meta_data_Set_Tests.tc,Meta_data_Set_Tests.Shape);
%%
path2colormap = strcat('Utilities\ScientificColourMaps8\','lipari','\','lipari','.mat');

load(path2colormap);

cmap        = colormap(lipari);

path2colormap = strcat('Utilities\ScientificColourMaps8\','glasgow','\','glasgow','.mat');

load(path2colormap);

cmap2        = colormap(glasgow);


% Figure
fg1 = scatter_plot_post_process();
clf;
z =(Data_S.Lambda0);
a = Data_S.tc./365.25./24./60./60./1e6./Data_S.n;
b = ((Data_S.tdet)).*Data_S.tc./Data_S.n;
b = b./365.25./24./60./60./1e6;
fg1.figure_number = 2;
fg1.x = a;
fg1.y = b;
fg1.c = log10(z);
fg1.xlim =[10^-1,10^2];
fg1.ylim =[10^-1,10^4];
%fg1.colormap_discrete = 20;
fg1.logx = 'log';
fg1.logy = 'log';
fg1.logcolor = 0;
fg1.colormap_f = 'lipari';
fg1.xlabel = '$t_c$';
fg1.ylabel = '$t_0$';
fg1.clim = [-4,0];
%fg2.clabel = '${\nu_{M,0}}$';
fg1.save_path = ptsave;% ptsave;
fg1.name_figure = 'tanal_t0d';
fg1.size_picture = [12.5,13];
fg1.make_scatter_plot;


allowed=a>10^(-1) & a<2*10^2;
%%
figure(2)
fg2 = scatter_plot_post_process();
colormap_f = [34, 139, 34;240, 230, 140;218, 165, 32;220, 20, 60]./255;
clf;
z =1./Data_S.tdet;
y = (Data_S.tau_max.*Data_S.tau0)./1e9;
x = Data_S.Lambda0;
z = Data_S.L0;
x = x(allowed);
y = y(allowed);
z = z(allowed);

fg2.figure_number = 2;
fg2.x = x;
fg2.y = y;
fg2.c = z;
fg2.xlim =[10^-4,10];
fg2.colormap_discrete = 10;
fg2.logx = 'log';
fg2.logy = 'linear';
fg2.logcolor = 0;
fg2.loop_unique_c = 1;
fg2.colormap_f = colormap_f;
fg2.xlabel = '$\Lambda_0$';
fg2.ylabel = '$\tau [GPa]$';
%fg2.clabel = '${\nu_{M,0}}$';
fg2.save_path = ptsave;% ptsave;
fg2.name_figure = 'Lambda_0_vs_tau_dim';
fg2.size_picture = [12.5,13];
fg2.make_scatter_plot;

%%
fg3 = scatter_plot_post_process();

z =Data_S.tdet;
y = (Data_S.tau_max);
x = Data_S.Lambda0;
z = Data_S.xdisl;
x = x(allowed);
y = y(allowed);
z = z(allowed);
fg3.figure_number = 3;
fg3.x = x;
fg3.y = y;
fg3.c = [0.0,0.0,0.0];
fg3.logx = 'log';
fg3.logy = 'linear';
fg3.logcolor = 0;
fg3.xlim =[10^-4,10];
fg3.size = 10;
%fg3.colormap_f = [0,0,0];
fg3.xlabel = '$\Lambda_0$';
fg3.ylabel = '$\tau^{\dagger}$';
fg3.save_path = ptsave;% ptsave;
fg3.name_figure = 'Lambda_0_vs_tau_dnd';
fg3.size_picture = [12.5,13];
fg3.make_scatter_plot;

%%
% Prepare the 3D visualisation
time_d_dim  = td_test.*tc_test./Meta_data_Set_Tests.secMyear;
allowed     = suc_test(:) == 1 & time_d_dim(:)>0.1 & time_d_dim(:) <120.0;
not_        = allowed==0 & ~isnan(Lambda_test(:));


figure(3)
clf;
ax = gca;
p1=scatter3(Vn_test(allowed),Vd_test(allowed),TS_test(allowed)-273.15,30,log10(Lambda_test(allowed)),'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
hold on
p2 = scatter3(Vn_test(not_),Vd_test(not_),TS_test(not_)-273.15,6,log10(Lambda_test(not_)),'filled','d','MarkerFaceAlpha',0.2,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
colormap(lipari)

ax.LineWidth = 1.2;
ax.Box = 'on';
view([-45 60 30])
caxis([-4,0]);
colorbar;
fname = ['3D_TS_Plot.png'];
pt = fullfile(ptsave,fname);
print(fname,'-dpng')
figure(3)
colormap(ax,lipari);
fname = ['CBar_1.png'];
pt = fullfile(fname);
print(fname,'-dpng')
%%

figure(4)
clf;
ax = gca;
p1=scatter3(Vn_test(allowed),Vd_test(allowed),L0_tests(allowed)./1e3,30,log10(Lambda_test(allowed)),'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.7);
hold on
p2 = scatter3(Vn_test(not_),Vd_test(not_),L0_tests(not_)./1e3,6,log10(Lambda_test(not_)),'filled','d','MarkerFaceAlpha',0.2,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
colormap(lipari);

ax.LineWidth = 1.2;
ax.Box = 'on';
%ax.BoxStyle = 'full';
view([-45 60 30])
caxis([-4,0]);
colorbar;
fname = ['3D_L0_Plot.png'];
pt = fullfile(fname);
print(fname,'-dpng')
%%
figure(5)
clf;
ax = gca;
p1=scatter(Vn_test(allowed),TS_test(allowed)-273.15,30,log10(Lambda_test(allowed)),'filled','MarkerEdgeColor','k');
hold on
p2 = scatter(Vn_test(not_),TS_test(not_)-273.15,5,log10(Lambda_test(not_)),'filled','d','MarkerFaceAlpha',0.1,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
colormap(lipari);
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
xline(15e-6)
xline(20e-6)
ax.LineWidth = 1.2;
caxis([-5,-1]);
colorbar;
fname = ['2D_TS_Plot.png'];
pt = fullfile(fname);
print(fname,'-dpng')

allowed2 = suc_test(:) == 1 & time_d_dim(:)>0.1 & time_d_dim(:) <120.0 & Vn_test(:) >=15e-6 & Vn_test(:) <=20e-6;




