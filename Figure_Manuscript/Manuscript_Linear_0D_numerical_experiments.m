%=========================================================================%
% Figure Linear 0D numerical experiments
%======================================================================
close all;
clear all;
clf;
addpath Utilities
tic;
load('../../Data_Base/Linear_Tests_Data_Base.mat','n_3','Data_S');
toc;

%% 
font_axes = 16; 
font_legend = 14; 
font_text   = 5; 
size_shit = [12,12.5];
LineWidth = 1.0; 
marker_size = 10;
folder = '../../Manuscript/'
figure_folder = 'figure_2';
supplementary_folder = 'figure_2_SUP'; 
folder_save=fullfile(folder,figure_folder);
folder_supplementary=fullfile(folder,supplementary_folder);

if not(isdir(folder_save))
    mkdir(folder_save);
end

if not(isdir(folder_supplementary))
    mkdir(folder_supplementary);
end
% Figure 1
% The first figure depicts the dimensional evolution of D, with dimensional
% time expressed in Myrs. And the corresponding detachment timescales
% nearby plotted against log10 Lambda.Then the same counterpart plots in
% dimensionless units.
%%
fun_0D = Manuscript_function_Container; 

[Fig1_A]=fun_0D.select_tests_prepare_variables(n_3,50.0,'time','D_norm','none','Linear',1,'xius');
[Fig1_C]=fun_0D.select_tests_prepare_variables(n_3,50.0,'time_nd','dDdt2','Lambda','Linear',1,'xius');

%% Figure Karlsruhe Presentation 



F3 = figure(1);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

a = squeeze(Fig1_C(1,:,:));
b = squeeze(Fig1_C(2,:,:));
c = squeeze(Fig1_C(3,:,:));

F1C.size_tests = length(squeeze(Fig1_C(1,1,:)));
F1C.color_lists = fun_0D.color_computation(F1C.size_tests,c,-7,0);
F1C.cmap = colormap(crameri('bilbao'));
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
F1C.c = colorbar(gca,'southoutside');
F1C.x.TickLabelInterpreter = 'latex';
F1C.y.TickLabelInterpreter = 'latex';
F1C.z.TickLabelInterpreter = 'latex';
F1C.c.Label.String = '$log_{10}\left(\Lambda\right)$';
F1C.c.Label.Interpreter = 'latex';
F1C.c.TickLabelInterpreter = 'latex';
F1C.c.FontSize = font_axes; 
% Define coloraxis
F1C.c.Limits =[F1C.c_lim(1), F1C.c_lim(2)];
% min/max are rounded using the log values of Lambda
Ticks=round((F1C.c_lim(1))):round((F1C.c_lim(2)));
% Produce the tick
TickLabels=arrayfun(@(x) sprintf('$10^{%d}$',x),Ticks,'UniformOutput',false);
%TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting

caxis([-7 0])
grid on
grid minor
box on
view([-30,20,45])

filename = 'FigureA';

pt=fullfile(folder_save,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')

F5=figure(2)
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
a = Data_S.Lambda;
b = Data_S.tdet;
c = Data_S.xiUS;
d = Data_S.n;
color_marker = Data_S.tc./(365.25*24*60*60*1e6);
a = a(c==50.0 & d==3.5);
b = b(c==50.0 & d==3.5);
color_marker = color_marker(c==50.0 & d==3.5);
hold on 
line([10^-6,1.0],[1.0,1.0],'Color',[.85 .85 .85],'LineWidth',1.0)
F1B.p = scatter(a,b,15,log10(color_marker),'filled');
colormap(crameri('berlin'));
F1B.c    = colorbar(gca,'southoutside');
F1B.c.Label.String = '$log_{10}\left(t_c\right) \mathrm{[Myrs]}$';
F1B.c.Label.Interpreter = 'latex';
F1B.c.Limits = [-3,2];
Ticks =  [-3,-2,-1,0,1,2];
Tickslabel = {'$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'};
F1B.c.Ticks = Ticks;
F1B.c.TickLabelInterpreter = 'latex';
F1B.c.TickLabels = Tickslabel;
F1B.c.FontSize = font_axes; 
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
F1B.x.Limits   = [10^-6,1.0];
%F1B.y.Limits   = [0.7,50];
F1B.x.FontSize=font_axes;
F1B.y.FontSize=font_axes;
F1B.x.TickLabelInterpreter = 'latex';
F1B.y.TickLabelInterpreter = 'latex';
F1B.x.Label.String = '$log_{10}\left(\Lambda\right)$';
F1B.x.Label.Interpreter = 'latex';
F1B.y.Label.String = '$t_d^{\dagger}$';
F1B.y.Label.Interpreter = 'latex';
F1B.x.MinorTick = 'on';
F1B.y.MinorTick = 'on';
set(gca,'Box','on','Layer','top')
grid on
grid minor
F1B.x.Scale = 'log';


filename = 'FigureB';

pt=fullfile(folder_save,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')
bla = 0;

%% Fig Supplementary Linear Manuscript
%==========================================================================
% Figure showing the effective stress evolution, and how t_d is correlated
% with the effective stress. The purpouse of the picture is showing how the
% stress dictactes the evolution of the experiments.
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

F3 = figure(4);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
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

filename = 'FigureSB1';

pt=fullfile(folder_supplementary,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')

F4 = figure(7);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
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
filename = 'FigureSB2';

pt=fullfile(folder_supplementary,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')

F5 = figure(5);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
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
filename = 'FigureSB3';
pt=fullfile(folder_supplementary,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')
F6 = figure(6);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

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
filename = 'FigureSB4';

pt=fullfile(folder_supplementary,filename);
set(gcf,'Color','w')
print(pt,'-r600','-dpng')



