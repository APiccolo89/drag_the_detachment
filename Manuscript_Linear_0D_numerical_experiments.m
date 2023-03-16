%=========================================================================%
% Figure Linear 0D numerical experiments
%======================================================================
close all;
clear all;
clf;
addpath Utilities
addpath \Users\'Andrea Piccolo'\Dropbox\freezeColors-master\
addpath \Users\Andrea Piccolo\Dropbox\Bayreuth_Marcel\export_fig
tic;
load('Linear_Tests_Data_Base.mat','n_3','Data_S');
toc;

%% 
font_axes = 8; 
font_legend = 8; 
font_text   = 5; 
size_shit = [13,13.5];
LineWidth = 0.8; 
marker_size = 6;

% Figure 1
% The first figure depicts the dimensional evolution of D, with dimensional
% time expressed in Myrs. And the corresponding detachment timescales
% nearby plotted against log10 Lambda.Then the same counterpart plots in
% dimensionless units.
%%
fun_0D = Manuscript_function_Container; 

[Fig1_A]=fun_0D.select_tests_prepare_variables(n_3,50.0,'time','D_norm','none','Linear',1,'xius');
[Fig1_C]=fun_0D.select_tests_prepare_variables(n_3,50.0,'time_nd','D_norm','Lambda','Linear',1,'xius');

%% Figure Karlsruhe Presentation 
figure(100)
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20,20], 'PaperUnits', 'centimeters', 'PaperSize', [20, 20])

a = Data_S.Lambda;
b = Data_S.tau_det;
b2 = Data_S.tau_max;

c = Data_S.xiUS;
d = Data_S.n;
color_marker = Data_S.tc./(365.25*24*60*60*1e6);
a = a(c==50.0 & d==3.5);
b = b(c==50.0 & d==3.5);
b2 = b2(c==50.0 & d==3.5);

color_marker = color_marker(c==50.0 & d==3.5);
hold on
F1B.p(1) = scatter(a,b,15,[0.8 0.8 0.8],'filled');
F1B.p(2) = scatter(a,b2,15,[153 204 0]./255,'filled');
l=legend('$\tau^{\dagger}_d$','$\tau^{\dagger}_{max}$');
l.Interpreter = 'latex'
line([10^(-6) 10^(0)],[1.0,1.0],'LineWidth',0.8,'Color',[0.8 0.8 0.8])
%colormap(crameri('devon'));freezeColors;
%F1B.p.MarkerFaceColor = 'k';
F1B.x = get(gca,'XAxis');
F1B.y = get(gca,'YAxis');
F1B.x.Color = [0.0 0.0 0.0];
F1B.x.LineWidth = LineWidth;
F1B.y.Color = [0.0 0.0 0.0];
F1B.y.LineWidth = LineWidth;
F1B.x.Limits   = [10^-6,1];
F1B.y.Limits   = [0.5,20];
F1B.x.FontSize=16;
F1B.y.FontSize=16;
F1B.x.TickLabelInterpreter = 'latex';
F1B.y.TickLabelInterpreter = 'latex';
F1B.x.Label.String = '$log_{10}\left(\Lambda\right)$';
F1B.x.Label.Interpreter = 'latex';
F1B.y.Label.String = '$\tau_{eff}^{\dagger} $';
F1B.y.Label.Interpreter = 'latex';
F1B.x.MinorTick = 'on';
F1B.y.MinorTick = 'on';
set(gca,'Box','on','Layer','top')
grid on
grid minor
F1B.x.Scale = 'log';
F1B.y.Scale = 'linear';
filename = 'Figure_Linear_1_scatterTau_ND';
folder   = 'D:\Dropbox\Bayreuth_Marcel\Conference_Abstract_Poster_Presentation\Presentation\Poster_KARLSRUHE_PIC';

if not(isdir(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')


%%
F=figure(1);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
nexttile
a = squeeze(Fig1_A(1,:,:));
b = squeeze(Fig1_A(2,:,:));
F1A.p = plot(a,b,'Color','#8467c1','LineWidth',LineWidth);
F1A.x = get(gca,'XAxis');
F1A.y = get(gca,'YAxis');
F1A.x.Color = [0.0 0.0 0.0];
F1A.x.LineWidth = LineWidth;
F1A.y.Color = [0.0 0.0 0.0];
F1A.y.LineWidth = LineWidth;
F1A.x.Limits   = [0,350];
F1A.y.Limits   = [0.1,1.01];
F1A.x.FontSize=font_axes;
F1A.y.FontSize=font_axes;
F1A.x.TickLabelInterpreter = 'latex';
F1A.y.TickLabelInterpreter = 'latex';
F1A.x.Label.String = '$t \mathrm{[Myrs]}$';
F1A.x.Label.Interpreter = 'latex';
F1A.y.Label.String = '$D^{\dagger}$';
F1A.y.Label.Interpreter = 'latex';
F1A.x.MinorTick = 'on';
F1A.y.MinorTick = 'on';
grid on
grid minor
set(gca,'Box','on','Layer','top')
%a = annotation(1.01,0.99,'a','FontSize',font_text,'Units','normalized',Interpreter='latex',FontWeight='bold',FontName='helvetica',)
F1A.a = annotation('textbox','LineWidth',1.2);
F1A.pos_ = get(gca,'Position');
F1A.a.Position = [F1A.pos_(1)+F1A.pos_(3)+0.01,F1A.pos_(2)+0.99*F1A.pos_(4),0.001,0.0005]; 
F1A.a.String = 'a'; 
F1A.a.FontSize = font_text;
F1A.a.BackgroundColor = 'k'; 
F1A.a.FaceAlpha = 0.3; 
F1A.a.FitBoxToText = 'on';
F1A.a.FontWeight = 'bold'; 

F1A.a.Color = 'w'; 

nexttile;
a = Data_S.Lambda;
b = Data_S.tc.*Data_S.tdet./(365.25*24*60*60*1e6);
c = Data_S.xiUS;
d = Data_S.n;
color_marker = Data_S.tc./(365.25*24*60*60*1e6);
a = a(c==50.0 & d==3.5);
b = b(c==50.0 & d==3.5);
color_marker = color_marker(c==50.0 & d==3.5);
F1B.p = scatter(a,b,15,log10(color_marker),'filled');
colormap(crameri('devon'));freezeColors;
F1B.c    = colorbar(gca,'southoutside');
F1B.c.Label.String = '$log_{10}\left(t_c\right) \mathrm{[Myrs]}$';
F1B.c.Label.Interpreter = 'latex';
F1B.c.Limits = [-3,2];
Ticks =  [-3,-2,-1,0,1,2];
Tickslabel = {'$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'};
F1B.c.Ticks = Ticks;
F1B.c.TickLabelInterpreter = 'latex';
F1B.c.TickLabels = Tickslabel;
F1B.c.FontSize = font_text; 
F1B.c.Color    = [0,0,0];
F1B.p.SizeData = 10;
F1B.p.Marker = 'o';
%F1B.p.MarkerFaceColor = 'k';
F1B.x = get(gca,'XAxis');
F1B.y = get(gca,'YAxis');
F1B.x.Color = [0.0 0.0 0.0];
F1B.x.LineWidth = LineWidth;
F1B.y.Color = [0.0 0.0 0.0];
F1B.y.LineWidth = LineWidth;
F1B.x.Limits   = [10^-6,0.1];
F1B.y.Limits   = [0.1,430];
F1B.x.FontSize=font_axes;
F1B.y.FontSize=font_axes;
F1B.x.TickLabelInterpreter = 'latex';
F1B.y.TickLabelInterpreter = 'latex';
F1B.x.Label.String = '$log_{10}\left(\Lambda\right)$';
F1B.x.Label.Interpreter = 'latex';
F1B.y.Label.String = '$t_d \mathrm{[Myrs]}$';
F1B.y.Label.Interpreter = 'latex';
F1B.x.MinorTick = 'on';
F1B.y.MinorTick = 'on';
set(gca,'Box','on','Layer','top')
grid on
grid minor
F1B.x.Scale = 'log';
F1B.y.Scale = 'log';

freezeColors; freezeColors(F1B.c);
F1B.a = annotation('textbox','LineWidth',1.2);
F1B.pos_ = get(gca,'Position');
F1B.a.Position = [F1B.pos_(1)+F1B.pos_(3)+0.01,F1B.pos_(2)+0.99*F1B.pos_(4),0.001,0.0005]; 
F1B.a.String = 'b'; 
F1B.a.FontSize = font_text;
F1B.a.BackgroundColor = 'k'; 
F1B.a.FaceAlpha = 0.3; 
F1B.a.FitBoxToText = 'on';
F1B.a.FontWeight = 'bold'; 

F1B.a.Color = 'w'; 
nexttile
a = squeeze(Fig1_C(1,:,:));
b = squeeze(Fig1_C(2,:,:));
c = squeeze(Fig1_C(3,:,:));

F1C.size_tests = length(squeeze(Fig1_C(1,1,:)));
F1C.color_lists = fun_0D.color_computation(F1C.size_tests,c,-7,0);
F1C.cmap = colormap(crameri('bilbao'));freezeColors;
F1C.c_lim = [-7,0];
for i = 1:length(c(1,:))
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    F1C.p(i)=plot3(aa,bb,cc,'Color',F1C.cmap(F1C.color_lists(i),:),'LineWidth',LineWidth);
    F1C.p2(i)=plot3(aa,bb,ones(size(aa)).*10^-10,'Color','k','LineWidth',LineWidth,'LineStyle','-');
end
hold off

F1C.x = get(gca,'XAxis');
F1C.y = get(gca,'YAxis');
F1C.z = get(gca,'ZAxis');
F1C.x.FontSize = font_axes;
F1C.y.FontSize = font_axes;
F1C.z.FontSize = font_axes;
F1C.x.Label.String = '$t^{\dagger}$';
F1C.x.Label.Interpreter = 'latex';
F1C.y.Label.String = '$D^{\dagger}$';
F1C.y.Label.Interpreter = 'latex';
F1C.z.Label.String = '$\Lambda$';
F1C.z.Label.Interpreter = 'latex';
F1C.x.Color = [0.0 0.0 0.0];
F1C.x.LineWidth = LineWidth;
F1C.y.Color = [0.0 0.0 0.0];
F1C.y.LineWidth = LineWidth;
F1C.z.Color = [0.0 0.0 0.0];
F1C.z.LineWidth = 1.0;
F1C.x.MinorTick = 'on';
F1C.y.MinorTick = 'on';
F1C.z.MinorTick = 'on';
F1C.z.Scale = 'log';
F1C.c = colorbar;
F1C.x.TickLabelInterpreter = 'latex';
F1C.y.TickLabelInterpreter = 'latex';
F1C.z.TickLabelInterpreter = 'latex';
F1C.c.Label.String = '$log_{10}\left(\Lambda\right)$';
F1C.c.Label.Interpreter = 'latex';
F1C.c.TickLabelInterpreter = 'latex';
% Define coloraxis
F1C.c.Limits =[F1C.c_lim(1), F1C.c_lim(2)];
% min/max are rounded using the log values of Lambda
Ticks=round((F1C.c_lim(1))):round((F1C.c_lim(2)));
% Produce the tick
TickLabels=arrayfun(@(x) sprintf('$10^{%d}$',x),Ticks,'UniformOutput',false);
%TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
F1C.c.Ticks= Ticks;
F1C.c.TickLabelInterpreter = 'latex';
F1C.c.TickLabels=TickLabels;
F1C.c.Location = 'southoutside';
F1C.c.FontSize = font_text; 
F1C.c.Color    = [0,0,0];
freezeColors(F1C.c)
caxis([-7 0])
grid on
grid minor
box on
view([-30,20,45])
F1C.a = annotation('textbox','LineWidth',1.2);
F1C.pos_ = get(gca,'Position');
F1C.a.Position = [F1C.pos_(1)+F1C.pos_(3)+0.01,F1C.pos_(2)+0.99*F1C.pos_(4),0.001,0.0005]; 
F1C.a.String = 'c'; 
F1C.a.FontSize = font_text;
F1C.a.BackgroundColor = 'k'; 
F1C.a.FaceAlpha = 0.3; 
F1C.a.FitBoxToText = 'on';
F1C.a.FontWeight = 'bold'; 

F1C.a.Color = 'w'; 
nexttile;

a = Data_S.Lambda;
b = Data_S.tdet;
c = Data_S.xiUS;
d = Data_S.n;
a = a(c==50.0 & d==3.5);
b = b(c==50.0 & d==3.5);
F1D.p = scatter(a,b,marker_size,'filled');
F1D.p.MarkerFaceColor = '#8467c1';
F1D.x = get(gca,'XAxis');
F1D.y = get(gca,'YAxis');
F1D.x.Color = [0.0 0.0 0.0];
F1D.x.LineWidth = 1.2;
F1D.y.Color = [0.0 0.0 0.0];
F1D.y.LineWidth = 1.2;
F1D.x.Limits   = [10^-6,0.1];
%F1D.y.Limits   = [0.1,430];
F1D.x.FontSize=font_axes;
F1D.y.FontSize=font_axes;
F1D.x.TickLabelInterpreter = 'latex';
F1D.y.TickLabelInterpreter = 'latex';
F1D.x.Label.String = '$log_{10}\left(\Lambda\right)$';
F1D.x.Label.Interpreter = 'latex';
F1D.y.Label.String = '$t_d^{\dagger}$';
F1D.y.Label.Interpreter = 'latex';
F1D.x.MinorTick = 'on';
F1D.y.MinorTick = 'on';
set(gca,'Box','on','Layer','top')
grid on
grid minor
F1D.x.Scale = 'log';
F1D.y.Scale = 'linear';
F1D.y.Limits = [0.8,60.0];
F1D.x.Limits = [10^(-7),1];
F1D.a = annotation('textbox','LineWidth',1.2);
F1D.pos_ = get(gca,'Position');
F1D.a.Position = [F1D.pos_(1)+F1D.pos_(3)+0.01,F1D.pos_(2)+0.99*F1D.pos_(4),0.001,0.0005]; 
F1D.a.String = 'd'; 
F1D.a.FontSize = font_text;
F1D.a.BackgroundColor = 'k'; 
F1D.a.FaceAlpha = 0.3; 
F1D.a.FitBoxToText = 'on';
F1D.a.FontWeight = 'bold'; 

F1D.a.Color = 'w'; 
filename = 'Figure_1';
folder   = '../Manuscript_Main';
if not(isfolder(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf, 'Color', 'w') 
export_fig(pt,'-r600')
bla = 0;

%% Fig Supplementary Linear Manuscript
%==========================================================================
% Figure showing the effective stress evolution, and how t_d is correlated
% with the effective stress. The purpouse of the picture is showing how the
% stress dictactes the evolution of the experiments.
%==========================================================================
F2 = figure(2);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

[Fig2_A]=fun_0D.select_tests_prepare_variables(n_3,50.0,'time_nd','tau_eff','none','Linear',1);
[Fig2_B]=fun_0D.select_tests_prepare_variables(n_3,50.0,'tau_B','tauD','Lambda','Linear',1);
[Fig2_C]=fun_0D.select_tests_prepare_variables(n_3,50.0,'time_nd','Drag','Lambda','Linear',1);

a = squeeze(Fig2_A(1,:,:));
b = squeeze(Fig2_A(2,:,:));
tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
nexttile
F2A.p = plot(a,b,'Color','#e4ba1e','LineStyle','-','LineWidth',LineWidth);
F2A.x = get(gca,'XAxis');
F2A.y = get(gca,'YAxis');
F2A.x.Color = [0.0 0.0 0.0];
F2A.x.LineWidth = LineWidth;
F2A.y.Color = [0.0 0.0 0.0];
F2A.y.LineWidth = LineWidth;
F2A.x.Limits   = [0,50];
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
F2A.a = annotation('textbox','LineWidth',1.2);
F2A.pos_ = get(gca,'Position');
F2A.a.Position = [F2A.pos_(1)+F2A.pos_(3)+0.01,F2A.pos_(2)+0.99*F2A.pos_(4),0.001,0.0005]; 
F2A.a.String = 'a'; 
F2A.a.FontSize = font_text;
F2A.a.BackgroundColor = 'k'; 
F2A.a.FaceAlpha = 0.3; 
F2A.a.FitBoxToText = 'on';
F2A.a.FontWeight = 'bold'; 
F2A.a.Color     = 'w'; 

nexttile
a = squeeze(Fig2_B(1,:,:));
b = squeeze(Fig2_B(2,:,:));
c = squeeze(Fig2_B(3,:,:));
F2B.color_lists = fun_0D.color_computation(length(c(1,:)),c,-7,0);
F2B.cmap = colormap(crameri('bilbao'));freezeColors;
colormap(crameri('bilbao'));freezeColors;
F2B.c_lim = [-7,0];
hold on
for i=1:length(c(1,:))
    F2B.p(i)=plot(a(:,i),b(:,i),'Color',F2B.cmap(F2B.color_lists(i),:),'LineWidth',LineWidth);
end
F2B.x = get(gca,'XAxis');
F2B.y = get(gca,'YAxis');
F2B.x.Color = [0.0 0.0 0.0];
F2B.x.LineWidth = LineWidth;
F2B.y.Color = [0.0 0.0 0.0];
F2B.y.LineWidth = LineWidth;
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
F2B.c.Label.String = '$log_{10}\left(\Lambda\right)$';
F2B.c.Label.Interpreter = 'latex';
F2B.c.TickLabelInterpreter = 'latex';
% Define coloraxis
F2B.c.Limits =[F2B.c_lim(1), F2B.c_lim(2)];
% min/max are rounded using the log values of Lambda
Ticks=round((F2B.c_lim(1))):round((F2B.c_lim(2)));
% Produce the tick
TickLabels=arrayfun(@(x) sprintf('$10^{%d}$',x),Ticks,'UniformOutput',false);
%TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
F2B.c.Ticks= Ticks;
F2B.c.TickLabelInterpreter = 'latex';
F2B.c.TickLabels=TickLabels;
F2B.c.Location = 'southoutside';
F2B.c.FontSize = font_text; 
F2B.c.Color = [0 0 0 ]; 
caxis([-7,0]);
freezeColors(F2B.c)
F2B.a = annotation('textbox','LineWidth',1.2);
F2B.pos_ = get(gca,'Position');
F2B.a.Position = [F2B.pos_(1)+F2B.pos_(3)+0.01,F2B.pos_(2)+0.99*F2B.pos_(4),0.001,0.0005]; 
F2B.a.String = 'b'; 
F2B.a.FontSize = font_text;
F2B.a.BackgroundColor = 'k'; 
F2B.a.FaceAlpha = 0.3; 
F2B.a.FitBoxToText = 'on';
F2B.a.FontWeight = 'bold'; 
F2B.a.Color     = 'w'; 
nexttile;
a = Data_S.Lambda;
b = Data_S.tdet;
c = Data_S.xiUS;
d = Data_S.n;
e = Data_S.tau_det;
f = Data_S.tau_max;
a = a(c==50.0 & d==3.5);
b = b(c==50.0 & d==3.5);
e = e(c==50.0 & d==3.5);
f = f(c==50.0 & d==3.5);
%x1 axis
hold on
F2C.p1 = scatter(a,f,marker_size);
F2C.p1.Marker = 'diamond';
F2C.p1.MarkerFaceColor = '#749d82';
F2C.p1.MarkerEdgeColor = 'k';
F2C.p2 = scatter(a,e,marker_size);
F2C.p2.Marker = 'o';
F2C.p2.MarkerFaceColor = '#4d6ba2';
F2C.p2.MarkerEdgeColor = 'k';
hold off
F2C.x = get(gca,'XAxis');
F2C.y = get(gca,'YAxis');
F2C.x.Color = [0.0 0.0 0.0];
F2C.x.LineWidth = LineWidth;
F2C.y.Color = [0.0 0.0 0.0];
F2C.y.LineWidth = LineWidth;
%F2B.x.Limits   = [0,50];
%F2B.y.Limits   = [0.4,20.00];
F2C.x.FontSize=font_axes;
F2C.y.FontSize=font_axes;
F2C.x.TickLabelInterpreter = 'latex';
F2C.y.TickLabelInterpreter = 'latex';
F2C.x.Label.String = '$\Lambda$';
F2C.x.Label.Interpreter = 'latex';
F2C.y.Label.String = '$\tau^{\dagger}$';
F2C.y.Label.Interpreter = 'latex';
F2C.x.Scale = 'log';
F2C.x.MinorTick = 'on';
F2C.y.MinorTick = 'on';
F2C.l = legend('$\tau^{\dagger \mathrm{MAX}}_{eff}$','$\tau^{\dagger \mathrm{det}}_{eff}$');
F2C.l.Interpreter = 'latex';
F2C.l.FontSize = font_legend;
grid on
grid minor
set(gca,'Box','on','Layer','top')
F2C.a = annotation('textbox','LineWidth',1.2);
F2C.pos_ = get(gca,'Position');
F2C.a.Position = [F2C.pos_(1)+F2C.pos_(3)+0.01,F2C.pos_(2)+0.99*F2C.pos_(4),0.001,0.0005]; 
F2C.a.String = 'c'; 
F2C.a.FontSize = font_text;
F2C.a.BackgroundColor = 'k'; 
F2C.a.FaceAlpha = 0.3; 
F2C.a.FitBoxToText = 'on';
F2C.a.FontWeight = 'bold'; 
F2C.a.Color     = 'w'; 

nexttile;
a = squeeze(Fig2_C(1,:,:));
b = squeeze(Fig2_C(2,:,:));
c = squeeze(Fig2_C(3,:,:));
F2D.color_lists = fun_0D.color_computation(length(c(1,:)),c,-7,0);
F2D.cmap = colormap(crameri('bilbao'));freezeColors;
colormap(crameri('bilbao'));freezeColors;
F2D.c_lim = [-7,0];
hold on
for i=1:length(c(1,:))
    F2D.p(i)=plot(a(:,i),1-b(:,i),'Color',F2B.cmap(F2B.color_lists(i),:),'LineWidth',LineWidth);
end
F2D.x = get(gca,'XAxis');
F2D.y = get(gca,'YAxis');
F2D.x.Color = [0.0 0.0 0.0];
F2D.x.LineWidth = LineWidth;
F2D.y.Color = [0.0 0.0 0.0];
F2D.y.LineWidth = LineWidth;
%FD.x.Limits   = [0,50];
%F2B.y.Limits   = [0.4,20.00];
F2D.x.FontSize=font_axes;
F2D.y.FontSize=font_axes;
F2D.x.TickLabelInterpreter = 'latex';
F2D.y.TickLabelInterpreter = 'latex';
F2D.x.Label.String = '$t^{\dagger}$';
F2D.x.Label.Interpreter = 'latex';
F2D.y.Label.String = '$F_{eff}^{\dagger}$';
F2D.y.Label.Interpreter = 'latex';
F2D.x.Scale = 'linear';
F2D.y.Scale = 'linear';
F2D.x.MinorTick = 'on';
F2D.y.MinorTick = 'on';
grid on
grid minor

set(gca,'Box','on','Layer','top')
F2D.a = annotation('textbox','LineWidth',1.2);
F2D.pos_ = get(gca,'Position');
F2D.a.Position = [F2D.pos_(1)+F2D.pos_(3)+0.01,F2D.pos_(2)+0.99*F2D.pos_(4),0.001,0.0005]; 
F2D.a.String = 'd'; 
F2D.a.FontSize = font_text;
F2D.a.BackgroundColor = 'k'; 
F2D.a.FaceAlpha = 0.3; 
F2D.a.FitBoxToText = 'on';
F2D.a.FontWeight = 'bold';
F2D.a.Color     = 'w'; 

caxis([-7,0]);
filename = 'Stress_Linear_Supplementary';
folder   = '../Manuscript_Supplementary';
if not(isfolder(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
 set(gcf, 'Color', 'w') 
export_fig(pt,'-r600')

bla = 0;
%% Second figure of Linear Mantle
% Figure showing the functional relation between n, xius with t_d
%=========================================================================
x = Data_S.Lambda;
y = Data_S.tdet;
n = Data_S.n;
xius = Data_S.xiUS;
unique_n = unique(n);
unique_xius = unique(xius);
c_n = 3.5;
c_xi=50.0;

F3 = figure(3);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

tiledlayout(2,2,'Padding','compact','TileSpacing','compact')
nexttile
hold on
color_code = {'#8cafc2','#b8c28c','#ad7f7f','#f0cba8'};
marker_type     = {'d','o','square','*'};
for i =1: length(unique_n)
    n_p = unique_n(i);
    ch = xius == c_xi & n == n_p;
    x_p = x(ch==1);
    y_p = y(ch==1);
    F3A.p(i) = scatter(x_p,y_p,marker_size,'filled');
    F3A.p(i).Marker = marker_type{i};
    F3A.p(i).MarkerFaceColor = color_code{i};
    F3A.p(i).MarkerEdgeColor ='k';
    F3A.p(i).LineWidth  = 0.5;
    string_legend{i} = (['$n = ',num2str(n_p),'$']);
end
hold off
F3A.l = legend(string_legend);
F3A.l.Interpreter = 'latex';
F3A.l.FontSize = font_legend;
F3A.l.Location = 'northwest';
F3A.x = get(gca,'XAxis');
F3A.y = get(gca,'YAxis');
F3A.x.LineWidth = LineWidth;
F3A.x.Scale = 'log';
F3A.x.Color = [0,0,0];
%F3A.x.TickLabelInterpreter = 'latex';
%F3A.x.Label.String = '$\Lambda$';
%F3A.x.Label.Interpreter = 'latex';
F3A.x.TickLabels =[]; 
F3A.y.LineWidth = LineWidth;
F3A.y.Color = [0,0,0];
F3A.y.TickLabelInterpreter = 'latex';
F3A.y.Label.String = '$t_d^{\dagger}$';
F3A.y.Label.Interpreter = 'latex';
set(gca,'Box','on','Layer','top')
F3A.x.FontSize = font_axes;
F3A.y.FontSize = font_axes;

grid on
grid minor
F3A.a = annotation('textbox','LineWidth',LineWidth);
F3A.pos_ = get(gca,'Position');
F3A.a.Position = [F3A.pos_(1)+F3A.pos_(3)+0.001,F3A.pos_(2)+0.99*F3A.pos_(4),0.001,0.0005]; 
F3A.a.String = 'a'; 
F3A.a.FontSize = font_text;
F3A.a.BackgroundColor = 'k'; 
F3A.a.FaceAlpha = 0.3; 
F3A.a.FitBoxToText = 'on';
F3A.a.FontWeight = 'bold';
F3A.a.Color     = 'w'; 
nexttile
hold on
color_code = {'#8cafc2','#b8c28c','#ad7f7f','#f0cba8'};
marker_type     = {'d','o','square','*'};
for i =1: length(unique_xius)
    xiup = unique_xius(i);
    ch = xius == xiup & n == c_n;
    x_p = x(ch==1);
    y_p = y(ch==1);
    F3B.p(i) = scatter(x_p,y_p,marker_size,'filled');
    F3B.p(i).Marker = marker_type{i};
    F3B.p(i).MarkerFaceColor = color_code{i};
    F3B.p(i).MarkerEdgeColor ='k';
    F3B.p(i).LineWidth  = 0.5;
    string_legend{i} = (['$\xi^{\mathrm{S}} = ',num2str(xiup),'$']);
end
hold off
F3B.l = legend(string_legend);
F3B.l.Interpreter = 'latex';
F3B.l.Location = 'northwest'; 
F3B.l.FontSize = font_legend;

F3B.x = get(gca,'XAxis');
F3B.y = get(gca,'YAxis');
F3B.x.LineWidth = LineWidth;
F3B.x.Scale = 'log';
F3B.x.Color = [0,0,0];
%F3B.x.TickLabelInterpreter = 'latex';
%F3B.x.Label.String = '$\Lambda$';
%F3B.x.Label.Interpreter = 'latex';
F3B.x.TickLabels =[]; 
F3B.y.LineWidth = LineWidth;
F3B.y.Color = [0,0,0];
F3B.y.TickLabelInterpreter = 'latex';
F3B.y.Label.String = '$t_d^{\dagger}$';
F3B.y.Label.Interpreter = 'latex';
set(gca,'Box','on','Layer','top')
grid on
grid minor% Picture as a function of n @ constant xius = 50.0; 
F3B.x.FontSize = font_axes;
F3B.y.FontSize = font_axes;

F3B.a = annotation('textbox','LineWidth',1.2);
F3B.pos_ = get(gca,'Position');
F3B.a.Position = [F3B.pos_(1)+F3B.pos_(3)+0.005,F3B.pos_(2)+0.99*F3B.pos_(4),0.001,0.0005]; 
F3B.a.String = 'b'; 
F3B.a.FontSize = font_text;
F3B.a.BackgroundColor = 'k'; 
F3B.a.FaceAlpha = 0.3; 
F3B.a.FitBoxToText = 'on';
F3B.a.FontWeight = 'bold';
F3B.a.Color     = 'w'; 

nexttile 
hold on
color_code = {'#8cafc2','#b8c28c','#ad7f7f','#f0cba8'};
marker_type     = {'d','o','square','*'};
for i =1: length(unique_n)
    n_p = unique_n(i);
    ch = xius == c_xi & n == n_p;
    x_p = x(ch==1);
    y_p = 1./y(ch==1);
    F3C.p(i) = scatter(x_p,y_p,marker_size,'filled');
    F3C.p(i).Marker = marker_type{i};
    F3C.p(i).MarkerFaceColor = color_code{i};
    F3C.p(i).MarkerEdgeColor ='k';
    F3C.p(i).LineWidth  = 0.5;
    string_legend{i} = (['$n = ',num2str(n_p),'$']);
end
hold off
F3C.l = legend(string_legend);
F3C.l.Interpreter = 'latex';
F3C.l.FontSize = font_legend;
F3C.l.Location = 'southwest';
F3C.x = get(gca,'XAxis');
F3C.y = get(gca,'YAxis');
F3C.x.LineWidth = LineWidth;
F3C.x.Scale = 'log';
F3C.x.Color = [0,0,0];
F3C.x.TickLabelInterpreter = 'latex';
F3C.x.Label.String = '$\Lambda$';
F3C.x.Label.Interpreter = 'latex';
%F3B.x.TickLabels =[]; 
F3C.y.LineWidth = LineWidth;
F3C.y.Color = [0,0,0];
F3C.y.TickLabelInterpreter = 'latex';
F3C.y.Label.String = '$\frac{1}{t_d^{\dagger}}$';
F3C.y.Label.Interpreter = 'latex';
F3C.y.Limits = [-0.01 1.1]; 
F3C.x.FontSize = font_axes;
F3C.y.FontSize = font_axes;

set(gca,'Box','on','Layer','top')
grid on
grid minor
F3C.a = annotation('textbox','LineWidth',1.2);
F3C.pos_ = get(gca,'Position');
F3C.a.Position = [F3C.pos_(1)+F3C.pos_(3)+0.001,F3C.pos_(2)+0.99*F3C.pos_(4),0.001,0.0005]; 
F3C.a.String = 'c'; 
F3C.a.FontSize = font_text;
F3C.a.BackgroundColor = 'k'; 
F3C.a.FaceAlpha = 0.3; 
F3C.a.FitBoxToText = 'on';
F3C.a.FontWeight = 'bold';
F3C.a.Color     = 'w'; 
nexttile
hold on
color_code = {'#8cafc2','#b8c28c','#ad7f7f','#f0cba8'};
marker_type     = {'d','o','square','*'};
for i =1: length(unique_xius)
    xiup = unique_xius(i);
    ch = xius == xiup & n == c_n;
    x_p = x(ch==1);
    y_p = 1./y(ch==1);
    F3D.p(i) = scatter(x_p,y_p,marker_size,'filled');
    F3D.p(i).Marker = marker_type{i};
    F3D.p(i).MarkerFaceColor = color_code{i};
    F3D.p(i).MarkerEdgeColor ='k';
    F3D.p(i).LineWidth  = 0.5;
    string_legend{i} = (['$\xi^{\mathrm{S}} = ',num2str(xiup),'$']);
end
hold off
F3D.l = legend(string_legend);
F3D.l.Interpreter = 'latex';
F3D.l.Location = 'southwest'; 
F3D.l.FontSize = font_legend;
F3D.x = get(gca,'XAxis');
F3D.y = get(gca,'YAxis');
F3D.x.FontSize = font_axes;
F3D.y.FontSize = font_axes;

F3D.x.LineWidth = LineWidth;
F3D.x.Scale = 'log';
F3D.x.Color = [0,0,0];
F3D.x.TickLabelInterpreter = 'latex';
F3D.x.Label.String = '$\Lambda$';
F3D.x.Label.Interpreter = 'latex';
%F3B.x.TickLabels =[]; 
F3D.y.LineWidth = LineWidth;
F3D.y.Color = [0,0,0];
F3D.y.TickLabelInterpreter = 'latex';
F3D.y.Label.String = '$\frac{1}{t_d^{\dagger}}$';
F3D.y.Label.Interpreter = 'latex';
F3D.y.Limits = [-0.01 1.1]; 
set(gca,'Box','on','Layer','top')
grid on
grid minor% Picture as a function of n @ constant xius = 50.0; 
F3D.a = annotation('textbox','LineWidth',1.2);
F3D.pos_ = get(gca,'Position');
F3D.a.Position = [F3D.pos_(1)+F3D.pos_(3)+0.005,F3D.pos_(2)+0.99*F3D.pos_(4),0.001,0.0005]; 
F3D.a.String = 'd'; 
F3D.a.FontSize = font_text;
F3D.a.BackgroundColor = 'k'; 
F3D.a.FaceAlpha = 0.3; 
F3D.a.FitBoxToText = 'on';
F3D.a.FontWeight = 'bold';
F3D.a.Color     = 'w'; 
bla = 0.0; 
filename = 'Effects_xi_n';
pt=fullfile(folder,filename);
 set(gcf, 'Color', 'w') 
export_fig(pt,'-r600')



