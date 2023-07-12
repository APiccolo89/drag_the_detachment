%=========================================================================%
%  Manuscript Non Linear Experiments 
%=========================================================================%
%

%=========================================================================%
% Figure Linear 0D numerical experiments
%======================================================================
close all;
clear all;
clf;
% Add relevant Path
addpath ../Realistic_Case/
addpath ../Utilities/
% load the data base from the folder of data base
load('../NonLinear_Tests_Data_Base_n3.mat','n_3','Data_S');

%%
% General information about the picture: s
font_axes = 16; 
font_legend = 14; 
font_text   = 5; 
size_picture = [12,12.5];
LineWidth = 1.0; 
marker_size = 10;
folder   = '../../Manuscript/';
% Compute before hand the t_diff to filter the experiment 
t_diff = 80e3^2./1e-6;
% call the classes that handle the Data structures to have the data that
% should be plot
fun_0D = Manuscript_function_Container;  

%% 
% Figure 1
% Notes: 
% a) The slab viscosity contrast has little effects for values higher than
% 100. 
% b) The biggest effects are concerned with the viscosity contrast of the
% upper mantle (i.e. the viscosity contrast between 

%% Figure Main 3
% Function that creates the data structure with the relevant information
% for the picture at hand 
[Fig1_A]=fun_0D.select_tests_prepare_variables(n_3,0,'time_nd','tau_eff','Lambda','NonLinear',3,'xius');
[Fig1_B]=fun_0D.select_tests_prepare_variables(n_3,0,'D_norm','tau_eff','Lambda','NonLinear',3,'xius');
[Fig1_C]=fun_0D.select_tests_prepare_variables(n_3,0,'dDdt','tauD_B','Lambda','NonLinear',3,'xius');


%%
cmap        = colormap(crameri('berlin')); 
% folder picture
% figure_folder = 'figure_3';
% supplementary_folder = 'figure_3_SUP'; 
% folder_save=fullfile(folder,figure_folder);
% folder_supplementary=fullfile(folder,supplementary_folder);
 size_tests = length(squeeze(Fig1_B(1,1,:)));
% 
% % if not(isdir(folder_save))
%     mkdir(folder_save);
% end
% 
% if not(isdir(folder_supplementary))
%     mkdir(folder_supplementary);
% end
F=figure(1);
clf; 
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])
ax_2 = axes; 
% First part
ax_2.Box = 'on';
a = squeeze(Fig1_B(1,:,:));
b = squeeze(Fig1_B(2,:,:));
c = squeeze(Fig1_B(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c,-6,0);

for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = 1-aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p2(i)=plot(aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',LineWidth);
end
ax_2.TickLabelInterpreter = 'latex';
ax_2.LineWidth = LineWidth; 
ax_2.Box = 'on';
ax_2.FontSize = font_axes;
% X Axis
ax_2.XColor = [0,0,0];
ax_2.XMinorTick = 'on';
ax_2.XMinorGrid = 'on';
ax_2.XGrid    = 'on'; 
ax_2.XLabel.Interpreter = 'latex';
ax_2.XLabel.FontWeight = 'bold'; 
% Y Axis 
ax_2.YColor =[0,0,0]; 
ax_2.XColor = [0,0,0];
ax_2.YMinorTick = 'on';
ax_2.YMinorGrid = 'on';
ax_2.YGrid    = 'on'; 
ax_2.YLabel.Interpreter = 'latex';
ax_2.YLabel.FontWeight = 'bold'; 
ax_2.YScale = 'log'; 
% Labels,Scales
ax_2.XLabel.String = '$1-D^{\dagger}$'; 
ax_2.YLabel.String = '$\tau^{\dagger}_{eff}$'; 
c    = colorbar(gca,'southoutside');
c.Label.Interpreter = 'latex';
c.Label.String = '$log_{10}\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$';
c.Limits = [-6,0];
Ticks =  [-6:1:0];

Tickslabel = {'{-6}$','${-5}$','${-4}$','${-3}$','$-2$','${-1}$','$0$'};
c.Ticks = Ticks;
c.TickLabelInterpreter = 'latex';
c.TickLabels = Tickslabel;
c.FontSize = font_axes; 
c.Color    = [0,0,0];
caxis([-6,0]); 

filename = 'Figure_3_A';
if not(isdir(folder_save))
    mkdir(folder_save);
end
pt=fullfile(folder_save,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')
% 
%%
F2=figure(2);

clf;
cmap        = colormap(crameri('berlin')); 

set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])
ax_2 = axes; 
% First part
ax_2.Box = 'on';
a = squeeze(Fig1_A(1,:,:));
b = squeeze(Fig1_A(2,:,:));
c = squeeze(Fig1_A(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c,-5,0);

for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p2(i)=plot(aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',LineWidth);
end
ax_2.TickLabelInterpreter = 'latex';
ax_2.LineWidth = LineWidth; 
ax_2.Box = 'on';
ax_2.FontSize = font_axes;
ax_2.XTickLabel = [];
ax_2.YTickLabel = [];

% X Axis
ax_2.XColor = [0,0,0];
ax_2.XMinorTick = 'on';
ax_2.XMinorGrid = 'on';
ax_2.XGrid    = 'on'; 
ax_2.XLabel.Interpreter = 'latex';
ax_2.XLabel.FontWeight = 'bold'; 
% Y Axis 
ax_2.YColor =[0,0,0]; 
ax_2.XColor = [0,0,0];
ax_2.YMinorTick = 'on';
ax_2.YMinorGrid = 'on';
ax_2.YGrid    = 'on'; 
ax_2.YLabel.Interpreter = 'latex';
ax_2.YLabel.FontWeight = 'bold'; 
ax_2.LineWidth = 1.2; 
ax_2.XColor = [0 0 0 ];
ax_2.YColor = [0 0 0 ];
ax_2.XMinorGrid = 'on';
ax_2.YMinorGrid = 'on';

% Labels,Scales
%ax_2.XLabel.String = '$t^{\dagger}$'; 
%ax_2.YLabel.String = '$\tau^{\dagger}_{\mathrm{eff}}$'; 
caxis([-5,0]); 

filename = 'Figure_3_StressSUP1';
if not(isdir(folder_supplementary))
    mkdir(folder_supplementary);
end
pt=fullfile(folder_supplementary,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')

delete(p2)
ax_2.Visible = 'off';
c    = colorbar(gca,'southoutside');
c.Label.Interpreter = 'latex';
%c.Label.String = '$log_{10}\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$';
c.Limits = [-5,0];
Ticks =  [-5:1:0];
%Tickslabel = {'${-5}$','${-4}$','${-3}$','$-2$','${-1}$','$0$'};
c.Ticks = Ticks;
c.TickLabelInterpreter = 'latex';
c.TickLabels = [];
c.FontSize = font_axes; 
c.Color    = [0,0,0];
ax_2.YScale = 'log'; 
filename = ['Figure_3_StressSUP1_CBAR'];

pt=fullfile(folder_supplementary,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')

%%

%%
F2=figure(2);
clf; 
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])

allowed = log10(Data_S.xiUM)>0; 
a = Data_S.Lambda(allowed==1)./(1+Data_S.xiUM(allowed==1)); 
b = 1./Data_S.tdet(allowed==1);
a2 =  Data_S.Lambda(allowed==0)./(1+Data_S.xiUM(allowed==0));
b2 = 1./Data_S.tdet(allowed==0);

ax_4 = axes;
hold on
p4 = scatter(a,b,marker_size,[30 88 143]./255,'filled','d','MarkerEdgeColor','k');
p5 = scatter(a2,b2,marker_size,[208 211 213]./255,'filled','o','MarkerEdgeColor','k');
ax_4.TickLabelInterpreter = 'latex';
ax_4.LineWidth = LineWidth; 
ax_4.Box = 'on';
ax_4.FontSize = font_axes;
% X4Axis
ax_4.XColor = [0,0,0];
ax_4.XMinorTick = 'on';
ax_4.XMinorGrid = 'on';
ax_4.XGrid    = 'on'; 
ax_4.XLabel.Interpreter = 'latex';
ax_4.XTickLabel =[];
% YAxis 
ax_4.YColor =[0,0,0]; 
ax_4.XColor = [0,0,0];
ax_4.YMinorTick = 'on';
ax_4.YMinorGrid = 'on';
ax_4.YGrid    = 'on'; 
ax_4.YLabel.Interpreter = 'latex';
ax_4.XLim = [10^-5,10^0];
ax_4.YLim = [0,1.1];
ax_4.YTickLabel =[];
ax_4.LineWidth = 1.2; 
ax_4.XColor = [0 0 0 ];
ax_4.YColor = [0 0 0 ];
ax_4.XMinorGrid = 'on';
ax_4.YMinorGrid = 'on';


% Labels,Scales
ax_4.XScale = 'log';
%ax_4.XLabel.String = '$\frac{\Lambda}{1+\xi^{\mathrm{UM}}}$'; 
%ax_4.YLabel.String = '$\frac{1}{t^{\dagger}_{\mathrm{d}}}$';

% Legend 
l = legend('$log_{10}\left(\xi^{UM}\right) > 0.0$','$log_{10}\left(\xi^{UM}\right) < 0.0$');
l.Interpreter = 'latex'; 
l.FontSize = font_legend;
l.Location = 'southwest';

filename = 'Figure_3_B';

pt=fullfile(folder_save,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')
%%
F3=figure(3);
clf; 
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])

allowed = log10(Data_S.xiUM)>0; 
a = Data_S.Lambda(allowed==1)./(1+Data_S.xiUM(allowed==1)); 
b = Data_S.tau_max(allowed==1);
a2 =  Data_S.Lambda(allowed==0)./(1+Data_S.xiUM(allowed==0));
b2 = Data_S.tau_max(allowed==0);

ax_4 = axes;
hold on
p4 = scatter(a,b,marker_size,[30 88 143]./255,'filled','d','MarkerEdgeColor','k');
p5 = scatter(a2,b2,marker_size,[208 211 213]./255,'filled','o','MarkerEdgeColor','k');
ax_4.TickLabelInterpreter = 'latex';
ax_4.LineWidth = LineWidth; 
ax_4.Box = 'on';
ax_4.FontSize = font_axes;
% X4Axis
ax_4.XColor = [0,0,0];
ax_4.XMinorTick = 'on';
ax_4.XMinorGrid = 'on';
ax_4.XGrid    = 'on'; 
ax_4.XLabel.Interpreter = 'latex';
% YAxis 
ax_4.XScale = 'log';

ax_4.YColor =[0,0,0]; 
ax_4.XColor = [0,0,0];
ax_4.YMinorTick = 'on';
ax_4.YMinorGrid = 'on';
ax_4.YGrid    = 'on'; 
ax_4.YLabel.Interpreter = 'latex';
ax_4.XLim = [10^-5,10^0];
ax_4.YLim = [0,10];
ax_4.XTickLabel =[];
ax_4.LineWidth = 1.2; 
ax_4.XColor = [0 0 0 ];
ax_4.YColor = [0 0 0 ];
ax_4.XMinorGrid = 'on';
ax_4.YMinorGrid = 'on';


% Labels,Scales
%ax_4.XLabel.String = '$\frac{\Lambda}{1+\xi^{\mathrm{UM}}}$'; 
%ax_4.YLabel.String = '$\frac{1}{t^{\dagger}_{\mathrm{d}}}$';
ax_4.YTickLabel =[];


% Legend 
l = legend('$log_{10}\left(\xi^{UM}\right) > 0.0$','$log_{10}\left(\xi^{UM}\right) < 0.0$');
l.Interpreter = 'latex'; 
l.FontSize = font_legend;
l.Location = 'southwest';

filename = 'Figure_3_C';

pt=fullfile(folder_save,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')

%% Figure Rheology for the supplementary
n=3.5;
xium = 10.^[-6:0.1:10]; 
tau0 = 100; 
tau = 0:1:2000; 
tau = tau./tau0; 
[X,t ]= meshgrid(xium,tau);

eta_eff = 1./(1+X.*t.^(n-1));
figure(3)
clf; 
ax = gca; 
p=pcolor(tau,xium,log10(eta_eff')); colormap(crameri('nuuk',10));shading interp;%colorbar;
ax.YScale = 'log';

hold on
xline(1.0,'Color','r',LineWidth=1.0)
yline(10,'LineStyle',':','LineWidth',1.2,'Color','k')
caxis([-10,0])
ax.XTickLabel = [];
ax.YTickLabel = []; 
ax.LineWidth = 1.2;
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
ax.Box ='on';
ax.Layer='top';
filename = 'Figure_Rheology_StressSUP1';
if not(isdir(folder_supplementary))
    mkdir(folder_supplementary);
end
pt=fullfile(folder_supplementary,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')

delete(p)
ax.Visible = 'off';
c    = colorbar(gca,'southoutside');
c.Label.Interpreter = 'latex';
%c.Label.String = '$log_{10}\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$';
c.Ticks = [];
c.TickLabelInterpreter = 'latex';
c.TickLabels = [];
c.FontSize = font_axes; 
c.Color    = [0,0,0];
filename = ['Rheology_CBar'];

pt=fullfile(folder_supplementary,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')




