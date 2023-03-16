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
load('NonLinear_Tests_Data_Base.mat','n_3','Data_S');
toc;
%%
font_axes = 8; 
font_legend = 8; 
font_text   = 8; 
size_shit = [13,14.5];
LineWidth = 0.8; 
marker_size = 6;

% Figure 1
% Notes: 
% a) The slab viscosity contrast has little effects for values higher than
% 100. 
% b) The biggest effects are concerned with the viscosity contrast of the
% upper mantle (i.e. the viscosity contrast between 
%%
t_diff = 80e3^2./1e-6;

fun_0D = Manuscript_function_Container; 
[Fig1_A]=fun_0D.select_tests_prepare_variables(n_3,t_diff,'time_nd','D_norm','xium','NonLinear',2,'xius');
%[Fig1_B]=fun_0D.select_tests_prepare_variables(n_3,t_diff,'time_nd','tau_eff','Lambda','NonLinear',2,'xius');
[Fig1_B]=fun_0D.select_tests_prepare_variables(n_3,t_diff,'time_nd','Psi','xium','NonLinear',2,'xius');
%%
size_tests = length(squeeze(Fig1_A(1,1,:)));
d_x = 5.18./size_shit(1);
d_y = 5.18./size_shit(2);
p_x = [0.8, 6.9 ]./size_shit(1)+0.05;
p_y    = [0.8, 7.0]./size_shit(2)+0.05;
F=figure(1);
clf; 
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
ax_1 = axes; 
a = squeeze(Fig1_A(1,:,:));
b = squeeze(Fig1_A(2,:,:));
c = squeeze(Fig1_A(3,:,:));

color_lists = fun_0D.color_computation(size_tests,c,2,10);
%cmap = colormap(crameri('berlin'));freezeColors;
c_lim = [2, 10];


for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p(i)=plot(ax_1,aa,bb,'Color','k','LineWidth',LineWidth);
end

ax_1.TickLabelInterpreter = 'latex';
ax_1.LineWidth = LineWidth; 
ax_1.Box = 'on';
ax_1.FontSize = font_text;
ax_1.Position=[p_x(1),p_y(2),d_x,d_y];
% X Axis
ax_1.XColor = [0,0,0];
ax_1.XMinorTick = 'on';
ax_1.XMinorGrid = 'on';
ax_1.XGrid    = 'on'; 
ax_1.XLabel.Interpreter = 'latex';
ax_1.XLabel.FontWeight = 'bold'; 
ax_1.YLim = [0.05,1.0];
% Y Axis 
ax_1.YColor =[0,0,0]; 
ax_1.XColor = [0,0,0];
ax_1.YMinorTick = 'on';
ax_1.YMinorGrid = 'on';
ax_1.YGrid    = 'on'; 
ax_1.YLabel.Interpreter = 'latex';
ax_1.YLabel.FontWeight = 'bold'; 
% Labels,Scales
ax_1.XLabel.String = '$t^{\dagger}$'; 
ax_1.YLabel.String = '$D^{\dagger}$'; 
% Second subplot
ax_2 = axes;
% First part
ax_2.Position = [p_x(2),p_y(2),d_x,d_y];
ax_2.Box = 'on';
a = squeeze(Fig1_B(1,:,:));
b = squeeze(Fig1_B(2,:,:));
c = squeeze(Fig1_B(3,:,:));
for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p2(i)=plot(ax_2,aa,bb,'Color','k','LineWidth',LineWidth);
end

ax_2.TickLabelInterpreter = 'latex';
ax_2.LineWidth = LineWidth; 
ax_2.Box = 'on';
ax_2.FontSize = font_text;
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
% Labels,Scales
ax_2.XLabel.String = '$t^{\dagger}$'; 
ax_2.YLabel.String = '$\frac{\Psi}{\Psi_0}$'; 
% Scatter Plot 
a = Data_S.Lambda./(1+Data_S.xiUM); 
b = Data_S.tdet;
c = Data_S.xiUM; 
ax_3 = axes;
ax_3.Position=[p_x(1),p_y(1),d_x,d_y];
p3 = scatter(ax_3,a,b,marker_size,'k','filled','d');

ax_3.CLim = [2,10]; 
ax_3.TickLabelInterpreter = 'latex';
ax_3.LineWidth = LineWidth; 
ax_3.Box = 'on';
ax_3.FontSize = font_text;
% X Axis
ax_3.XColor = [0,0,0];
ax_3.XMinorTick = 'on';
ax_3.XMinorGrid = 'on';
ax_3.XGrid    = 'on'; 
ax_3.XLabel.Interpreter = 'latex';
ax_3.XLabel.FontWeight = 'bold'; 
% Y Axis 
ax_3.YColor =[0,0,0]; 
ax_3.XColor = [0,0,0];
ax_3.YMinorTick = 'on';
ax_3.YMinorGrid = 'on';
ax_3.YGrid    = 'on'; 
ax_3.YLabel.Interpreter = 'latex';
ax_3.YLabel.FontWeight = 'bold'; 
% Labels,Scales
ax_3.XScale = 'log';
ax_3.XLabel.String = '$\frac{\Lambda}{1+\xi^{\mathrm{UM}}}$'; 
ax_3.YLabel.String = '$t^{\dagger}_d$'; 

a = Data_S.Lambda./(1+Data_S.xiUM); 
b = Data_S.tau_det;
c = Data_S.xiUM; 
ax_4 = axes;
ax_4.Position = [p_x(2),p_y(1),d_x,d_y];
p4 = scatter(ax_4,a,b,marker_size,'k','filled','d');

ax_4.CLim = [2,10]; 
ax_4.TickLabelInterpreter = 'latex';
ax_4.LineWidth = LineWidth; 
ax_4.Box = 'on';
ax_4.FontSize = font_text;
% X4Axis
ax_4.XColor = [0,0,0];
ax_4.XMinorTick = 'on';
ax_4.XMinorGrid = 'on';
ax_4.XGrid    = 'on'; 
ax_4.XLabel.Interpreter = 'latex';
ax_4.XLabel.FontWeight = 'bold'; 
% YAxis 
ax_4.YColor =[0,0,0]; 
ax_4.XColor = [0,0,0];
ax_4.YMinorTick = 'on';
ax_4.YMinorGrid = 'on';
ax_4.YGrid    = 'on'; 
ax_4.YLabel.Interpreter = 'latex';
ax_4.YLabel.FontWeight = 'bold'; 
% Labels,Scales
ax_4.XScale = 'log';
ax_4.XLabel.String = '$\frac{\Lambda}{1+\xi^{\mathrm{UM}}}$'; 
ax_4.YLabel.String = '$\tau^{\dagger det}_{\mathrm{eff}}$';
% Annotation
%%
a1 = annotation('textbox');
pos_ = ax_1.Position;
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
pos_ = ax_2.Position;
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
pos_ = ax_3.Position;
a3.Position = [pos_(1)+pos_(3)-0.055,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
a3.String = 'c'; 
a3.FontSize = font_text;
a3.BackgroundColor = 'k'; 
a3.FaceAlpha = 0.3; 
a3.LineWidth = 0.1;

a3.FitBoxToText = 'on';
a3.FontWeight = 'bold'; 
a3.Color     = 'w'; 


a4 = annotation('textbox');
pos_ = ax_4.Position;
a4.Position = [pos_(1)+pos_(3)-0.055,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
a4.String = 'd'; 
a4.FontSize = font_text;
a4.BackgroundColor = 'k'; 
a4.FaceAlpha = 0.3; 
a4.LineWidth = 0.1;
a4.FitBoxToText = 'on';
a4.FontWeight = 'bold'; 
a4.Color     = 'w'; 



bla = 0.0;

filename = 'Figure_Nonlinear_Psi_Psicrit';
folder   = '../Manuscript_Main';

if not(isfolder(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf,'Color','w')

export_fig(pt,'-r600')
bla = 0;

%% Figure Poster 
figure(5)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
ax_1 = axes; 
a = squeeze(Fig1_A(1,:,:));
b = squeeze(Fig1_A(2,:,:));
c = squeeze(Fig1_A(3,:,:));

color_lists = fun_0D.color_computation(size_tests,c,2,10);
cmap = colormap(crameri('berlin'));freezeColors;
c_lim = [2, 10];


for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p(i)=plot(ax_1,aa,bb,'Color',cmap(i,:),'LineWidth',LineWidth);
end

ax_1.TickLabelInterpreter = 'latex';
ax_1.LineWidth = LineWidth; 
ax_1.Box = 'on';
ax_1.FontSize = 16;
%ax_1.Position=[p_x(1),p_y(2),d_x,d_y];
% X Axis
ax_1.XColor = [0,0,0];
ax_1.XMinorTick = 'on';
ax_1.XMinorGrid = 'on';
ax_1.XGrid    = 'on'; 
ax_1.XLabel.Interpreter = 'latex';
ax_1.XLabel.FontWeight = 'bold'; 
ax_1.YLim = [0.05,1.0];
% Y Axis 
ax_1.YColor =[0,0,0]; 
ax_1.XColor = [0,0,0];
ax_1.YMinorTick = 'on';
ax_1.YMinorGrid = 'on';
ax_1.YGrid    = 'on'; 
ax_1.YLabel.Interpreter = 'latex';
ax_1.YLabel.FontWeight = 'bold'; 
% Labels,Scales
ax_1.XLabel.String = '$t^{\dagger}$'; 
ax_1.YLabel.String = '$D^{\dagger}$'; 
filename = 'Figure_Non_Linear_1_D';
folder   = 'D:\Dropbox\Bayreuth_Marcel\Conference_Abstract_Poster_Presentation\Presentation\Poster_KARLSRUHE_PIC';

if not(isdir(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')
bla = 0;
%%
figure(6)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
ax_1 = axes; 
a = squeeze(Fig1_B(1,:,:));
b = squeeze(Fig1_B(2,:,:));
c = squeeze(Fig1_B(3,:,:));

color_lists = fun_0D.color_computation(size_tests,c,-6,2);
cmap = colormap(parula(256));
c_lim = [-6, 2];


for i = 1:length(c(1,:))
    hold on 
    aa = a(:,i);
    bb = b(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    hold on
    p(i)=plot(ax_1,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',LineWidth);
end

ax_1.TickLabelInterpreter = 'latex';
ax_1.LineWidth = LineWidth; 
ax_1.Box = 'on';
ax_1.FontSize = 16;
%ax_1.Position=[p_x(1),p_y(2),d_x,d_y];
% X Axis
ax_1.XColor = [0,0,0];
ax_1.XMinorTick = 'on';
ax_1.XMinorGrid = 'on';
ax_1.XGrid    = 'on'; 
ax_1.XLabel.Interpreter = 'latex';
ax_1.XLabel.FontWeight = 'bold'; 
ax_1.YLim = [10^(-6),1.0];
ax_1.YScale = 'log';
% Y Axis 
ax_1.YColor =[0,0,0]; 
ax_1.XColor = [0,0,0];
ax_1.YMinorTick = 'on';
ax_1.YMinorGrid = 'on';
ax_1.YGrid    = 'on'; 
ax_1.YLabel.Interpreter = 'latex';
ax_1.YLabel.FontWeight = 'bold'; 
% Labels,Scales
ax_1.XLabel.String = '$t^{\dagger}$'; 
ax_1.YLabel.String = '$\Psi$'; 
c = colorbar; 
c.Label.String = '$log_{10}\left(\frac{\xi^{\mathrm{S}}}{\xi^{\mathrm{UM}}}\right)$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
% Define coloraxis
c.Limits =[-6, 2];
% min/max are rounded using the log values of Lambda
Ticks=round((-6:2));
% Produce the tick
TickLabels=arrayfun(@(x) sprintf('$10^{%d}$',x),Ticks,'UniformOutput',false);
%TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
c.Ticks= Ticks;
c.TickLabelInterpreter = 'latex';
c.TickLabels=TickLabels;
%c.Location = 'eastoutside';
c.FontSize = 12; 
c.Color    = [0,0,0];
caxis([-6,2])

filename = 'Figure_Non_Linear_1_Psi';
folder   = 'D:\Dropbox\Bayreuth_Marcel\Conference_Abstract_Poster_Presentation\Presentation\Poster_KARLSRUHE_PIC';

if not(isdir(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')
%%
figure(7)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
a = Data_S.Lambda./(1+Data_S.xiUM); 
b = Data_S.tau_max;
c = Data_S.xiUS; 

ax_3 = axes;
hold on
p3(1) = scatter(ax_3,a(Data_S.xiUS>10),b(Data_S.xiUS>10),marker_size,'k','filled','d');
p3(2) = scatter(ax_3,a(Data_S.xiUS<10),b(Data_S.xiUS<10),marker_size,'g','filled','d');
hold off
l = legend('$\xi^{S}> 10$','$\xi^{S} \leq 10$')
l.Interpreter = 'latex'
ax_3.TickLabelInterpreter = 'latex';
ax_3.LineWidth = LineWidth; 
ax_3.Box = 'on';
ax_3.FontSize = 16;
% X Axis
ax_3.XColor = [0,0,0];
ax_3.XMinorTick = 'on';
ax_3.XMinorGrid = 'on';
ax_3.XGrid    = 'on'; 
ax_3.XLabel.Interpreter = 'latex';
ax_3.XLabel.FontWeight = 'bold'; 
ax_3.XLim = [10^(-6),1];
% Y Axis 
ax_3.YColor =[0,0,0]; 
ax_3.XColor = [0,0,0];
ax_3.YMinorTick = 'on';
ax_3.YMinorGrid = 'on';
ax_3.YGrid    = 'on'; 
ax_3.YLabel.Interpreter = 'latex';
ax_3.YLabel.FontWeight = 'bold'; 
% Labels,Scales
ax_3.XScale = 'log';
ax_3.XLabel.String = '$\frac{\Lambda}{1+\xi^{\mathrm{UM}}}$'; 
ax_3.YLabel.String = '$\tau^{\dagger}_{max}$'; 
filename = 'Figure_Non_Linear_1_scatterTau';
folder   = 'D:\Dropbox\Bayreuth_Marcel\Conference_Abstract_Poster_Presentation\Presentation\Poster_KARLSRUHE_PIC';

if not(isdir(folder))
    mkdir(folder);
end
pt=fullfile(folder,filename);
set(gcf,'Color','w')

print(pt,'-r600','-dpng')

bla = 0;
%%

% nexttile(1,[3,1])
% a = squeeze(Fig1_A(1,:,:));
% b = squeeze(Fig1_A(2,:,:));
% F1A.gca= gca; 
% F1A.p = plot(a,b,'Color','#d1a282','LineWidth',LineWidth);
% F1A.x = get(gca,'XAxis');
% F1A.y = get(gca,'YAxis');
% F1A.x.Color = [0.0 0.0 0.0];
% F1A.x.LineWidth = LineWidth;
% F1A.y.Color = [0.0 0.0 0.0];
% F1A.x.LineWidth = LineWidth;
% F1A.x.Limits   = [0,15];
% %F1A.y.Limits   = [0.1,1.01];
% F1A.x.FontSize=font_axes;
% F1A.y.FontSize=font_axes;
% F1A.x.TickLabelInterpreter = 'latex';
% F1A.y.TickLabelInterpreter = 'latex';
% F1A.x.Label.String = '$t^{\dagger}$';
% F1A.x.Label.Interpreter = 'latex';
% F1A.y.Label.String = '$\frac{\Lambda}{\Lambda_{0}}$';
% F1A.y.Label.Interpreter = 'latex';
% F1A.x.MinorTick = 'on';
% F1A.y.MinorTick = 'on';
% F1A.y.Scale = 'log';
% grid on
% grid minor
% set(gca,'Box','on','Layer','top')
% %a = annotation(1.01,0.99,'a','FontSize',font_text,'Units','normalized',Interpreter='latex',FontWeight='bold',FontName='helvetica',)
% 
% 
% nexttile(2,[3,1])
% a = squeeze(Fig1_B(1,:,:));
% b = squeeze(Fig1_B(2,:,:));
% F1B.gca= gca; 
% F1B.p = plot(a,b,'Color','#1a9a3d','LineWidth',LineWidth);
% F1B.x = get(gca,'XAxis');
% F1B.y = get(gca,'YAxis');
% F1B.x.Color = [0.0 0.0 0.0];
% F1B.x.LineWidth = LineWidth;
% F1B.y.Color = [0.0 0.0 0.0];
% F1B.x.LineWidth = LineWidth;
% F1B.x.Limits   = [0,15];
% F1B.x.FontSize=font_axes;
% F1B.y.FontSize=font_axes;
% F1B.y.Scale = 'log';
% F1B.x.TickLabelInterpreter = 'latex';
% F1B.y.TickLabelInterpreter = 'latex';
% F1B.y.TickLabels = [];
% F1B.x.Label.String = '$t^{\dagger}$';
% F1B.x.Label.Interpreter = 'latex';
% F1B.y.Label.Interpreter = 'latex';
% F1B.x.MinorTick = 'on';
% F1B.y.MinorTick = 'on';
% F1B.y.Scale = 'log';
% 
% grid on
% grid minor
% set(gca,'Box','on','Layer','top')
% %a = annotation(1.01,0.99,'a','FontSize',font_text,'Units','normalized',Interpreter='latex',FontWeight='bold',FontName='helvetica',)
% 
% 
% nexttile(7,[2,1])
% a = Data_S.Lambda;
% b = Data_S.tdet;
% c = Data_S.xiUS;
% d = Data_S.xiUM; 
% a = a./(1+d); 
% d_u = unique(Data_S.xiUM); 
% hold on 
% color_code = {'#8cafc2','#b8c28c','#ad7f7f','#f0cba8','k','r'};
% F1C.gca= gca; 
% 
% for i = 1:length(d_u)
%     i = i 
%     d_up = d_u(i); 
%     chosen =c <=1e3 & d==d_up;
%     x = a(chosen==1);
%     y = 1./b(chosen==1);
% 
%     F1C.p(i) = scatter(x,y,marker_size,'filled','d');
%     F1C.p(i).MarkerFaceColor = 'k';
%     string_legend{i}= (['$log_{10}\left(\xi^{\mathrm{UM}}\right) = ',num2str(log10(d_up)),'$']);
%     chosen = []; 
% end
% 
% hold off
% F1C.x = get(gca,'XAxis');
% F1C.y = get(gca,'YAxis');
% F1C.x.Color = [0.0 0.0 0.0];
% F1C.x.LineWidth = LineWidth;
% F1C.y.Color = [0.0 0.0 0.0];
% F1C.x.LineWidth = LineWidth;
% F1C.x.Limits   = [10^-12,10^0];
% F1C.y.Limits   = [0.0,1.1];
% F1C.x.FontSize=font_axes;
% F1C.y.FontSize=font_axes;
% F1C.x.TickLabelInterpreter = 'latex';
% F1C.y.TickLabelInterpreter = 'latex';
% F1C.y.Label.String = '$\frac{1}{t_d^{\dagger}}$';
% F1C.y.Label.Interpreter = 'latex';
% F1C.x.MinorTick = 'on';
% F1C.y.MinorTick = 'on';
% F1C.x.Label.String = '$log_{10}\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$';
% F1C.x.Label.Interpreter = 'latex';
% 
% set(gca,'Box','on','Layer','top')
% grid on
% grid minor
% F1C.x.Scale = 'log';
% F1C.y.Scale = 'linear';
% 
% nexttile(8,[2,1])
% F1D.gca= gca; 
% 
% a = Data_S.Lambda;
% b = Data_S.tdet;
% c = Data_S.xiUS;
% d = Data_S.xiUM; 
% a = a./(1+d); 
% d_u = unique(Data_S.xiUM); 
% hold on 
% color_code = {'#8cafc2','#b8c28c','#ad7f7f','#f0cba8','k','r'};
% 
% for i = 1:length(d_u)
%     i = i 
%     d_up = d_u(i); 
%     chosen =c >=10^3 & d==d_up;
%     x = a(chosen==1);
%     y = 1./b(chosen==1);
% 
%     F1C.p(i) = scatter(x,y,marker_size,'filled','d');
%     F1C.p(i).MarkerFaceColor = 'k';
%     string_legend{i}= (['$log_{10}\left(\xi^{\mathrm{UM}}\right) = ',num2str(log10(d_up)),'$']);
%     chosen = []; 
% end
% 
% hold off
% F1D.x = get(gca,'XAxis');
% F1D.y = get(gca,'YAxis');
% F1D.x.Color = [0.0 0.0 0.0];
% F1D.x.LineWidth = LineWidth;
% F1D.y.Color = [0.0 0.0 0.0];
% F1D.x.LineWidth = LineWidth;
% F1D.x.Limits   = [10^-12,1.0];
% F1D.y.Limits   = [0.0,1.1];
% F1D.x.FontSize=font_axes;
% F1D.y.FontSize=font_axes;
% F1D.x.TickLabelInterpreter = 'latex';
% F1D.y.TickLabels = [];
% F1D.y.TickLabelInterpreter = 'latex';
% F1D.x.Label.String = '$log_{10}\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$';
% F1D.x.Label.Interpreter = 'latex';
% %F1D.y.Label.String = '$\frac{1}{t_d^{\dagger}}$';
% F1D.y.Label.Interpreter = 'latex';
% F1D.x.MinorTick = 'on';
% F1D.y.MinorTick = 'on';
% set(gca,'Box','on','Layer','top')
% grid on
% grid minor
% F1D.x.Scale = 'log';
% F1D.y.Scale = 'linear';
% 
% 
% % Annotation: 
% pos_=F1A.gca.Position;
% F1A.a = annotation('textbox','LineWidth',1.0);
% F1A.pos_ = get(gca,'Position');
% F1A.a.Position = [pos_(1)+pos_(3)+0.0001,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
% F1A.a.String = 'a'; 
% F1A.a.FontSize = font_text;
% F1A.a.BackgroundColor = 'k'; 
% F1A.a.FaceAlpha = 0.3; 
% F1A.a.FitBoxToText = 'on';
% F1A.a.FontWeight = 'bold'; 
% F1A.a.Color = 'w'; 
% 
% F1B.a = annotation('textbox','LineWidth',LineWidth);
% pos_ = F1B.gca.Position; 
% F1B.a.Position = [pos_(1)+pos_(3)+0.001,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
% F1B.a.String = 'b'; 
% F1B.a.FontSize = font_text;
% F1B.a.BackgroundColor = 'k'; 
% F1B.a.FaceAlpha = 0.3; 
% F1B.a.FitBoxToText = 'on';
% F1B.a.FontWeight = 'bold'; 
% F1B.a.Color = 'w';
% 
% F1C.a = annotation('textbox','LineWidth',LineWidth);
% pos_ = F1C.gca.Position;
% F1C.a.Position = [pos_(1)+pos_(3)+0.001,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
% F1C.a.String = 'c'; 
% F1C.a.FontSize = font_text;
% F1C.a.BackgroundColor = 'k'; 
% F1C.a.FaceAlpha = 0.3; 
% F1C.a.FitBoxToText = 'on';
% F1C.a.FontWeight = 'bold'; 
% F1C.a.Color = 'w'; 
% 
% F1D.a = annotation('textbox','LineWidth',LineWidth);
% pos_ = F1D.gca.Position;
% F1D.a.Position = [pos_(1)+pos_(3)+0.001,pos_(2)+0.99*pos_(4),0.001,0.0005]; 
% F1D.a.String = 'd'; 
% F1D.a.FontSize = font_text;
% F1D.a.BackgroundColor = 'k'; 
% F1D.a.FaceAlpha = 0.3; 
% F1D.a.FitBoxToText = 'on';
% F1D.a.FontWeight = 'bold'; 
% 
% F1D.a.Color = 'w'; 





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
