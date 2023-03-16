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
addpath Realistic_Case\ 
addpath Utilities
addpath \Users\'Andrea Piccolo'\Dropbox\freezeColors-master\
addpath \Users\Andrea Piccolo\Dropbox\Bayreuth_Marcel\export_fig
tic;
load('Realistic_DataSet.mat','Tests','Data_S');
toc;
%%
font_axes = 8; 
font_legend = 8; 
font_text   = 5; 
size_shit = [13,13.5];
LineWidth = 0.8; 
marker_size = 6;
n_3 = Tests; 
% Figure 1
% Notes: 
% a) The slab viscosity contrast has little effects for values higher than
% 100. 
% b) The biggest effects are concerned with the viscosity contrast of the
% upper mantle (i.e. the viscosity contrast between 
%%
fun_0D = Manuscript_function_Container; 
t_diff = 80e3^2./(1e-6);
[Fig1_A]=fun_0D.select_tests_prepare_variables(n_3,t_diff,'time_nd','D_norm','Lambda','NonLinear',2,'tcVdVn');
[Fig1_B]=fun_0D.select_tests_prepare_variables(n_3,t_diff,'time_nd','tau_eff','Lambda','NonLinear',2,'tcVdVn');
[Fig1_D]=fun_0D.select_tests_prepare_variables(n_3,t_diff,'time_nd','Psi','Lambda','NonLinear',2,'tcVdVn');





F = figure(1); 
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
tiledlayout(6,3,'Padding','compact','TileSpacing','compact')





%%
F=figure(1);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
tiledlayout(6,3,'Padding','compact','TileSpacing','compact')

nexttile(1,[5,1])
a = squeeze(Fig1_A(1,:,:));
b = squeeze(Fig1_A(2,:,:));
c = squeeze(Fig1_A(3,:,:));

size_tests = length(squeeze(Fig1_A(1,1,:)));
color_lists = fun_0D.color_computation(size_tests,c,-12,0);
cmap = colormap(crameri('bilbao'));freezeColors;
c_lim = [-12,0];


for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p(i)=plot(aa,bb,'Color','k','LineWidth',LineWidth);
end
ax = gca; 
ax.XLim = [0,40];
ax.YLim = [0.05,1.0];
ax.XMinorTick = 'on';
ax.XLabel.String = '$t^{\dagger}$';
ax.YLabel.String = '$D^{\dagger}$';
ax.YLabel.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.XColor = [0,0,0];
ax.YColor = [0,0,0]; 
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.LineWidth = 1.0; 
ax.Box = 'on';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = font_axes;



nexttile(2,[5,1])

a = squeeze(Fig1_B(1,:,:));
b = squeeze(Fig1_B(2,:,:));
c = squeeze(Fig1_B(3,:,:));

size_tests = length(squeeze(Fig1_A(1,1,:)));
color_lists = fun_0D.color_computation(size_tests,c,-12,0);
cmap = colormap(crameri('bilbao'));freezeColors;
c_lim = [-12,0];


for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p(i)=plot((aa),(bb),'Color','k','LineWidth',LineWidth);
end

ax2 = gca; 
ax2.XLim = [0,40];
ax2.YLim = [0.05,15];
ax2.XMinorTick = 'on';
ax2.XLabel.String = '$t^{\dagger}$';
ax2.YLabel.String = '$\tau_{eff}^{\dagger}$';
ax2.YLabel.Interpreter = 'latex';
ax2.XLabel.Interpreter = 'latex';
ax2.XColor = [0,0,0];
ax2.YColor = [0,0,0]; 
ax2.XGrid = 'on';
ax2.YGrid = 'on';
ax2.XMinorGrid = 'on';
ax2.YMinorGrid = 'on';
ax2.LineWidth = 1.0; 
ax2.Box = 'on';
ax2.TickLabelInterpreter = 'latex';
ax2.FontSize = font_axes;



nexttile(3,[5,1])

a = squeeze(Fig1_D(1,:,:));
b = squeeze(Fig1_D(2,:,:));
c = squeeze(Fig1_D(3,:,:));

size_tests = length(squeeze(Fig1_A(1,1,:)));
color_lists = fun_0D.color_computation(size_tests,c,-12,0);
cmap = colormap(crameri('bilbao'));freezeColors;


for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p(i)=plot((aa),(bb),'Color','k','LineWidth',LineWidth);
end
ax3 = (gca); 
ax3.XLim = [0,40];
%ax3.YLim = [0.05,1.0];
ax3.XMinorTick = 'on';
ax3.XLabel.String = '$t^{\dagger}$';
ax3.YLabel.String = '$\Psi$';
ax3.YLabel.Interpreter = 'latex';
ax3.XLabel.Interpreter = 'latex';
ax3.YScale = 'linear';
ax3.XColor = [0,0,0];
ax3.YColor = [0,0,0]; 
ax3.XGrid = 'on';
ax3.YGrid = 'on';
ax3.XMinorGrid = 'on';
ax3.YMinorGrid = 'on';
ax3.LineWidth = 1.0; 
ax3.Box = 'on';
ax3.TickLabelInterpreter = 'latex';
ax3.FontSize = font_axes;


a1 = annotation('textbox');
pos_ = ax.Position;
a1.Position = [pos_(1)+pos_(3)-0.055,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
a1.String = 'a'; 
a1.FontSize = font_text;
a1.BackgroundColor = 'k'; 
a1.FaceAlpha = 0.3; 
a1.LineWidth = 0.1;

a1.FitBoxToText = 'on';
a1.FontWeight = 'bold'; 
a1.Color     = 'w'; 

a2 = annotation('textbox');
pos_ = ax2.Position;
a2.Position = [pos_(1)+pos_(3)-0.055,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
a2.String = 'b'; 
a2.FontSize = font_text;
a2.BackgroundColor = 'k'; 
a2.FaceAlpha = 0.3; 
a2.LineWidth = 0.1;

a2.FitBoxToText = 'on';
a2.FontWeight = 'bold'; 
a2.Color     = 'w'; 


a3 = annotation('textbox');
pos_ = ax3.Position;
a3.Position = [pos_(1)+pos_(3)-0.055,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
a3.String = 'c'; 
a3.FontSize = font_text;
a3.BackgroundColor = 'k'; 
a3.FaceAlpha = 0.3; 
a3.LineWidth = 0.1;

a3.FitBoxToText = 'on';
a3.FontWeight = 'bold'; 
a3.Color     = 'w'; 



filename = 'Figure_Nonlinear_Realistic_Case2';
folder   = '../Manuscript_Main';

if not(isfolder(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf,'Color','w')

export_fig(pt,'-r600')


%%
%Figure 
% 
figure(2)
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

tc = Data_S.tc./Data_S.n;
ind = tc < t_diff; 
Lambda = Data_S.Lambda(ind==1)./(Data_S.xiUM(ind==1)+1); 
tdet   = Data_S.tdet(ind==1); 
t_slab = Data_S.T_Slab(ind==1); 
Psi    = Data_S.Psi(ind==1)./(1+Data_S.xiUM(ind==1)); 
taudet    = Data_S.tau_det(ind==1); 
taumax = Data_S.tau_max(ind==1); 
L0    = Data_S.L0(ind==1); 
Cd    = Data_S.Cd(ind==1); 
Cn    = Data_S.Cn(ind==1); 
phi   = Data_S.phi(ind==1); 
w     = Data_S.w(ind==1);
exp   = (w.*L0./(1+w.*L0.*phi)).*(-Cd+Cn);
Vd    = Data_S.Vd(ind==1);
Vn    = Data_S.Vn(ind==1);
eta_S_eff0 = Data_S.eta0DS./(1+Data_S.xiUS);
eta_UM_eff0 = Data_S.eta0DUM./(1+Data_S.xiUM); 
eta_S_eff0 = eta_S_eff0(ind==1);
eta_UM_eff0 = eta_UM_eff0(ind==1);
tiledlayout(4,2,'Padding','compact','TileSpacing','compact')

coord = [...
    2    15    0;
    10  15    0;
    10  20  0;
    2    20  0;
    2    15    1.0;
    10  15    1.0;
    10  20  1.0;
    2   20  1.0;];
idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; 3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';
nexttile(1,[3,1])
p = scatter3(Vd.*1e6,Vn.*1e6,1./tdet,20,(t_slab-273.15),'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.4); 
colormap(crameri('berlin',10));
xc = coord(:,1);
yc = coord(:,2);
zc = coord(:,3);
hold on 
l = patch(xc(idx), yc(idx), zc(idx), 'r', 'facealpha', 0.0,'LineWidth',1.0);
hold off 
view([-79 20])
ax = gca; 
ax.XLabel.String = '$V_\mathrm{n} [10^{-6}\frac{kJ}{m^3}]$';
ax.XLabel.Interpreter = 'latex'; 
ax.YLabel.String = '$V_\mathrm{n} [10^{-6}\frac{kJ}{m^3}]$';
ax.YLabel.Interpreter = 'latex'; 
ax.ZLabel.String = '$\frac{1}{t_d^{\dagger}}$';
ax.ZLabel.Interpreter = 'latex'; 
ax.TickLabelInterpreter ='latex'; 
ax.XMinorTick ='on';
ax.YMinorTick ='on';
ax.ZMinorTick ='on';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.ZMinorGrid = 'on';
ax.LineWidth = 1.0; 
ax.Box = 'on';

nexttile(2,[3,1])

coord = [...
    10^18    10^18    0;
    10^24  10^18    0;
    10^24  10^24  0;
    10^18   10^24  0;
    10^18    10^18    1.0;
    10^24  10^18    1.0;
    10^24  10^24  1.0;
    10^18  10^24  1.0;];

p = scatter3(eta_S_eff0,eta_UM_eff0,1./tdet,20,t_slab-273.15,'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.4); 
colormap(crameri('berlin',10))
xc = coord(:,1);
yc = coord(:,2);
zc = coord(:,3);
hold on 
l = patch(xc(idx), yc(idx), zc(idx), 'r', 'facealpha', 0.0,'LineWidth',1.0);
hold off 
view([-79 20])
ax = gca; 
ax.XLabel.String = '$\eta^{S}_{0,\mathrm{eff}}$';
ax.XLabel.Interpreter = 'latex'; 
ax.YLabel.String = '$\eta^{S}_{0,\mathrm{eff}}$';
ax.YLabel.Interpreter = 'latex'; 
ax.ZLabel.String = '$\frac{1}{t_d^{\dagger}}$';
ax.ZLabel.Interpreter = 'latex'; 
ax.TickLabelInterpreter ='latex'; 
ax.XScale = 'log';
ax.YScale = 'log';
ax.XMinorTick ='on';
ax.YMinorTick ='on';
ax.ZMinorTick ='on';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.ZMinorGrid = 'on';
ax.LineWidth = 1.0; 
ax.Box = 'on';
c = colorbar; 
c.Label.Interpreter = 'latex';
c.Label.String      = '$\bar{T}^{\mathrm{S}} ,[^{\circ}\mathrm{C}]$';
c.Color = [0,0,0];
c.Location = 'southoutside';
c.Position = [0.2,0.1,0.6,0.04];



filename = 'Figure_Nonlinear_Realistic_Case3';
folder   = '../Manuscript_Main';

if not(isfolder(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf,'Color','w')

export_fig(pt,'-transparent','-r600')
bla = 0;

%% Fig Supplementary Linear Manuscript
%==========================================================================
% Figure showing the effective stress evolution, and how t_d is correlated
% with the effective stress. The purpouse of the picture is showing how the
% stress dictactes the evolution of the experiments.
%==========================================================================
F2 = figure(2);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

[Fig2_A]=fun_0D.select_tests_prepare_variables(n_3,10^0,'time_nd','tau_eff','none','NonLinear',3);
[Fig2_B]=fun_0D.select_tests_prepare_variables(n_3,10^0,'tau_B','tauD','Lambda','NonLinear',3);
[Fig2_C]=fun_0D.select_tests_prepare_variables(n_3,10^0,'time_nd','Drag','Lambda','NonLinear',3);

a = squeeze(Fig2_A(1,:,:));
b = squeeze(Fig2_A(2,:,:));
tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
nexttile
F2A.gca = gca; 
F2A.p = plot(a,b,'Color','#e4ba1e','LineStyle','-','LineWidth',LineWidth);
F2A.x = get(gca,'XAxis');
F2A.y = get(gca,'YAxis');
F2A.x.Color = [0.0 0.0 0.0];
F2A.x.LineWidth = LineWidth;
F2A.y.Color = [0.0 0.0 0.0];
F2A.x.LineWidth = LineWidth;
F2A.x.Limits   = [0,25];
F2A.y.Limits   = [0.4,20.00];
F2A.x.FontSize=font_axes;
F2A.y.FontSize=font_axes;
F2A.x.TickLabelInterpreter = 'latex';
F2A.y.TickLabelInterpreter = 'latex';
F2A.x.Label.String = '$t^{\dagger}$';
F2A.x.Label.Interpreter = 'latex';
F2A.y.Label.String = '$\tau_{eff}^{\dagger}$';
F2A.y.Label.Interpreter = 'latex';
F2A.x.MinorTick = 'on';
F2A.y.MinorTick = 'on';
F2A.y.Scale = 'log';
grid on
grid minor
set(gca,'Box','on','Layer','top')

nexttile
a = squeeze(Fig2_B(1,:,:));
b = squeeze(Fig2_B(2,:,:));
c = squeeze(Fig2_B(3,:,:));
F2B.gca = gca; 

F2B.color_lists = fun_0D.color_computation(length(c(1,:)),c,-12,0);
F2B.cmap = colormap(crameri('bilbao'));freezeColors;
colormap(crameri('bilbao'));freezeColors;
F2B.c_lim = [-12,0];
hold on
for i=1:length(c(1,:))
    F2B.p(i)=plot(a(:,i),b(:,i),'Color',F2B.cmap(F2B.color_lists(i),:),'LineWidth',1.0);
end
F2B.x = get(gca,'XAxis');
F2B.y = get(gca,'YAxis');
F2B.x.Color = [0.0 0.0 0.0];
F2B.x.LineWidth = LineWidth;
F2B.y.Color = [0.0 0.0 0.0];
F2B.x.LineWidth = LineWidth;
%F2B.x.Limits   = [0,50];
%F2B.y.Limits   = [0.4,20.00];
F2B.x.FontSize=font_axes;
F2B.y.FontSize=font_axes;
F2B.x.TickLabelInterpreter = 'latex';
F2B.y.TickLabelInterpreter = 'latex';
F2B.x.Label.String = '$\tau_{B}^{\dagger}$';
F2B.x.Label.Interpreter = 'latex';
F2B.y.Label.String = '$\tau_{D}^{\dagger}$';
F2B.y.Label.Interpreter = 'latex';
F2B.x.Scale = 'log';
F2B.y.Scale = 'log';
F2B.x.MinorTick = 'on';
F2B.y.MinorTick = 'on';
grid on
grid minor
set(gca,'Box','on','Layer','top')
F2B.c = colorbar;
F2B.c.Label.String = '$log_{10}\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$';
F2B.c.Label.Interpreter = 'latex';
F2B.c.TickLabelInterpreter = 'latex';
% Define coloraxis
F2B.c.Limits =[F2B.c_lim(1), F2B.c_lim(2)];
% min/max are rounded using the log values of Lambda
Ticks=round((F2B.c_lim(1))):2:round((F2B.c_lim(2)));
% Produce the tick
TickLabels=arrayfun(@(x) sprintf('$10^{%d}$',x),Ticks,'UniformOutput',false);
%TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
F2B.c.Ticks= Ticks;
F2B.c.TickLabelInterpreter = 'latex';
F2B.c.TickLabels=TickLabels;
F2B.c.Location = 'southoutside';
F2B.c.FontSize = font_text; 
F2B.c.Color = [0 0 0 ]; 
caxis([-12,0]);
freezeColors(F2B.c)
 
nexttile;
a = Data_S.Lambda./(1+Data_S.xiUM);
b = Data_S.tdet;
c = Data_S.xiUS;
d = Data_S.n;
e = Data_S.tau_det;
f = Data_S.tau_max;
a = a(c>=1e1);
b = b(c>=1e1);
e = e(c>=1e1);
f = f(c>=1e1);
%x1 axis
hold on
F2C.gca = gca; 

F2C.p1 = scatter(a,f,20);
F2C.p1.Marker = 'diamond';
F2C.p1.MarkerFaceColor = '#749d82';
F2C.p1.MarkerEdgeColor = 'k';
F2C.p2 = scatter(a,e,15);
F2C.p2.Marker = 'o';
F2C.p2.MarkerFaceColor = '#4d6ba2';
F2C.p2.MarkerEdgeColor = 'k';
hold off
F2C.x = get(gca,'XAxis');
F2C.y = get(gca,'YAxis');
F2C.x.Color = [0.0 0.0 0.0];
F2C.x.LineWidth = LineWidth;
F2C.y.Color = [0.0 0.0 0.0];
F2C.x.LineWidth = LineWidth;
%F2B.x.Limits   = [0,50];
%F2B.y.Limits   = [0.4,20.00];
F2C.x.FontSize=font_axes;
F2C.y.FontSize=font_axes;
F2C.x.TickLabelInterpreter = 'latex';
F2C.y.TickLabelInterpreter = 'latex';
F2C.x.Label.String = '$\frac{\Lambda}{1+\xi^{\mathrm{UM}}}$';
F2C.x.Label.Interpreter = 'latex';
F2C.y.Label.String = '$\tau^{\dagger}$';
F2C.y.Label.Interpreter = 'latex';
F2C.x.Scale = 'log';
F2C.x.MinorTick = 'on';
F2C.y.MinorTick = 'on';
F2C.l = legend('$\tau^{\dagger \mathrm{MAX}}_{eff}$','$\tau^{\dagger \mathrm{det}}_{eff}$');
F2C.l.Interpreter = 'latex';
F2C.l.FontSize = font_legend;
F2C.l.Location = 'southwest';
grid on
grid minor
set(gca,'Box','on','Layer','top')

nexttile;
a = squeeze(Fig2_C(1,:,:));
b = squeeze(Fig2_C(2,:,:));
c = squeeze(Fig2_C(3,:,:));
F2D.gca = gca; 

F2D.color_lists = fun_0D.color_computation(length(c(1,:)),c,-12,0);
F2D.cmap = colormap(crameri('bilbao'));freezeColors;
colormap(crameri('bilbao'));freezeColors;
F2D.c_lim = [-12,0];
hold on
for i=1:length(c(1,:))
    F2D.p(i)=plot(a(:,i),1-b(:,i),'Color',F2D.cmap(F2D.color_lists(i),:),'LineWidth',1.0);
end
F2D.x = get(gca,'XAxis');
F2D.y = get(gca,'YAxis');
F2D.x.Color = [0.0 0.0 0.0];
F2D.x.LineWidth = LineWidth;
F2D.y.Color = [0.0 0.0 0.0];
F2D.x.LineWidth = LineWidth;
%FD.x.Limits   = [0,50];
%F2B.y.Limits   = [0.4,20.00];
F2D.x.FontSize=font_axes;
F2D.y.FontSize=font_axes;
F2D.x.TickLabelInterpreter = 'latex';
F2D.y.TickLabelInterpreter = 'latex';
F2D.x.Label.String = '$t^{\dagger}$';
F2D.x.Label.Interpreter = 'latex';
F2D.y.Label.String = '$F_{\mathrm{eff}}^{\dagger}$';
F2D.y.Label.Interpreter = 'latex';
F2D.x.Scale = 'linear';
F2D.y.Scale = 'log';
F2D.x.MinorTick = 'on';
F2D.y.MinorTick = 'on';
grid on
grid minor

set(gca,'Box','on','Layer','top')

caxis([-12,0]);

pos_=F2A.gca.Position;
F2A.a = annotation('textbox','LineWidth',1.0);
F2A.pos_ = get(gca,'Position');
F2A.a.Position = [pos_(1)+pos_(3)+0.0001,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
F2A.a.String = 'a'; 
F2A.a.FontSize = font_text;
F2A.a.BackgroundColor = 'k'; 
F2A.a.FaceAlpha = 0.3; 
F2A.a.FitBoxToText = 'on';
F2A.a.FontWeight = 'bold'; 
F2A.a.Color = 'w'; 

F2B.a = annotation('textbox','LineWidth',LineWidth);
pos_ = F2B.gca.Position; 
F2B.a.Position = [pos_(1)+pos_(3)+0.001,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
F2B.a.String = 'b'; 
F2B.a.FontSize = font_text;
F2B.a.BackgroundColor = 'k'; 
F2B.a.FaceAlpha = 0.3; 
F2B.a.FitBoxToText = 'on';
F2B.a.FontWeight = 'bold'; 
F2B.a.Color = 'w';

F2C.a = annotation('textbox','LineWidth',LineWidth);
pos_ = F2C.gca.Position;
F2C.a.Position = [pos_(1)+pos_(3)+0.001,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
F2C.a.String = 'c'; 
F2C.a.FontSize = font_text;
F2C.a.BackgroundColor = 'k'; 
F2C.a.FaceAlpha = 0.3; 
F2C.a.FitBoxToText = 'on';
F2C.a.FontWeight = 'bold'; 
F2C.a.Color = 'w'; 

F2D.a = annotation('textbox','LineWidth',LineWidth);
pos_ = F2D.gca.Position;
F2D.a.Position = [pos_(1)+pos_(3)+0.001,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
F2D.a.String = 'd'; 
F2D.a.FontSize = font_text;
F2D.a.BackgroundColor = 'k'; 
F2D.a.FaceAlpha = 0.3; 
F2D.a.FitBoxToText = 'on';
F2D.a.FontWeight = 'bold'; 

F2D.a.Color = 'w'; 




filename = 'Stress_Linear_Supplementary_Non_Linear';
folder   = '../Manuscript_Supplementary';
if not(isfolder(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf,'Color','w')
export_fig(pt,'-transparent','-r600')
%% 
clear all; 
close all; 
clf; 
