% Figure Manuscript Non Linear portion
clear all
close all
addpath('Utilities/')
addpath('Utilities/Plot_Class/')
addpath('Utilities/ScientificColourMaps8/')
clf; 
%%NonLinear_Tests_Data_Base_n3_iteration_high_T.mat It=load('NonLinear_Tests_Data_Base_n3_iteration.mat');/NonLinear_Tests_Data_Base_delta_L0
%It=load('..\Data_Base\NonLinear_Tests_Data_Base_n3_iteration_high_T_initial_guess.mat');

A = '../Data_Base/NonLinear_Regime_Main_data_base_gamma_1000e3.mat';
B = '../Data_Base/Linear_Regime_Main_data_base_gamma_1000e3.mat';
C = '../Data_Base/Transition_Regime_Main_data_base_gamma_1000e3.mat';
NL = load(A);
L = load(B);
T = load(C);
S1 = L.Data_S;
S2 = T.Data_S;
S3 = NL.Data_S;
fun_0D = Manuscript_function_Container;  

[Linear]=fun_0D.select_tests_prepare_variables(L.n_3,0,'time_nd','D_norm','Lambda0','NonLinear',3,'xius',[]);
[Transition]=fun_0D.select_tests_prepare_variables(T.n_3,0,'time_nd','D_norm','Lambda0','NonLinear',3,'xius',[]);
[Nonlinear]=fun_0D.select_tests_prepare_variables(NL.n_3,0,'time_nd','D_norm','Lambda0','NonLinear',3,'xius',[]);
path_save = 'New_Manuscript_Figure/';
%%  Lambda Figure 
close all;

f=figure(1);
clf;
set(gcf, 'Units','centimeters', 'Position', [0, 0, 15,15], 'PaperUnits', 'centimeters', 'PaperSize', [15, 15])
% linear experiment
subplot(2,3,1)
% linear experiment
ax0 = gca();
p0 = 0.08; 
ax0.Position = [0.08,0.65,0.28,0.3];

s1 = scatter(ax0,S1.Lambda0,1./S1.tdet,30,log10(S1.xiUM),"filled","MarkerEdgeColor",'k'); 
yline(1.0,'LineWidth',1.2,'Color','r','LineStyle',':')
[ax0] = set_default_Lambda(ax0);
ax0.YLabel.Position(1)=0.0000005;
ax0.YLabel.Position(2)=0.5;
caxis(ax0,[-12,12]);
colormap(ax0,crameri('oslo'));

subplot(2,3,2)
% lNinear experiment
ax1 = gca();
p2 = ax0.Position(1)+ax0.Position(3)+0.025;
ax1.Position = [p2,0.65,0.28,0.3];

s2 = scatter(ax1,S2.Lambda0,1./S2.tdet,30,log10(S2.xiUM),"filled","MarkerEdgeColor",'k'); 
yline(1.0,'LineWidth',1.2,'Color','r','LineStyle',':')
[ax1] = set_default_Lambda(ax1);
ax1.YTickLabel = [];
ax1.YLabel.String = [];
caxis(ax1,[-12,12]);
colormap(ax1,crameri('oslo'));

subplot(2,3,3)
% lNinear experiment
ax2 = gca();
p3 = ax1.Position(1)+ax1.Position(3)+0.025;

ax2.Position = [p3,0.65,0.28,0.3];
s3 = scatter(ax2,S3.Lambda0,1./S3.tdet,30,log10(S3.xiUM),"filled","MarkerEdgeColor",'k'); 
yline(1.0,'LineWidth',1.2,'Color','r','LineStyle',':')
[ax2] = set_default_Lambda(ax2);
ax2.YTickLabel = [];
ax2.YLabel.String = [];
caxis(ax2,[-12,12]);
colormap(ax2,crameri('oslo'));


% Function that creates the data structure with the relevant information
% for the picture at hand 
path2colormap = strcat('Utilities\ScientificColourMaps8\','lipari','\','lipari','.mat');

load(path2colormap);

cmap        = colormap(lipari);

size_tests = length(squeeze(Linear(1,1,:)));
c = squeeze(Linear(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c,-4,1);
subplot(2,3,4)
ax3 = gca(); 
ax3.Position = [p0,0.2,0.28,0.35];
x = squeeze(Linear(1,:,:));
y = squeeze(Linear(2,:,:));
c = squeeze(Linear(3,:,:));
for i = 1:size_tests
    hold on
    aa = x(:,i);
    bb = y(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    if ~isempty(aa)
        pl(i)=plot(ax3,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.9);
    end
    hold off
end
colormap(ax3,cmap);
ax3 = set_default_necking(ax3); 
ax3.YLabel.Position(1)=0.0000005;
ax3.YLabel.Position(2)=0.5;
caxis(ax3,[-4,1])




subplot(2,3,5)
ax4 = gca(); 
ax4.Position = [p2,0.2,0.28,0.35];


size_tests = length(squeeze(Transition(1,1,:)));
c = squeeze(Transition(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c,-4,1);

x = squeeze(Transition(1,:,:));
y = squeeze(Transition(2,:,:));
c = squeeze(Transition(3,:,:));



for i = 1:size_tests
    hold on
    aa = x(:,i);
    bb = y(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    if ~isempty(aa)
        pl2(i)=plot(ax4,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.9);
    end
    hold off
end

ax4 = set_default_necking(ax4); 

ax4.YTickLabel = [];
ax4.YLabel.String = [];
caxis(ax4,[-4,1])


subplot(2,3,6)
ax5 = gca(); 
ax5.Position = [p3,0.2,0.28,0.35];


size_tests = length(squeeze(Nonlinear(1,1,:)));
c = squeeze(Nonlinear(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c,-4,1);

x = squeeze(Nonlinear(1,:,:));
y = squeeze(Nonlinear(2,:,:));
c = squeeze(Nonlinear(3,:,:));



for i = 1:size_tests
    hold on
    aa = x(:,i);
    bb = y(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    if ~isempty(aa)
        pl3(i)=plot(ax5,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.9);
    end
    hold off
end

ax5 = set_default_necking(ax5); 
ax5.YTickLabel = [];
ax5.YLabel.String = [];
caxis(ax5,[-4,1])


% cbar = colorbar(ax2,"north");
% cbar.Box = 'on';
% cbar.Position = [0.22,0.85,0.60,0.02];
% cbar.Color = 'k';
% cbar.LineWidth = 1.2;
% cbar.TickLabelInterpreter = 'latex';
% cbar.Units = 'centimeters';
% cbar.Label.String = '$log_{10}\left(\xi_M\right)$';
% cbar.Label.Interpreter = 'latex';
% bla=0.0;

colormap(ax0,crameri('oslo'))
colormap(ax1,crameri('oslo'))
colormap(ax2,crameri('oslo'))

cbar=colorbar(ax0,"southoutside");
cbar.Position = [0.15,0.095,0.35,0.015];
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$log_{10}\left(\xi_M\right)$';
cbar.TickLabelInterpreter = 'latex';

cbar2= colorbar(ax5,"southoutside");
cbar2.Limits = [-4 1];
cbar2.Position = [0.55,0.095,0.35,0.015];

Ticks =  [-4,-3,-2,-1,0];
for i=1:numel(Ticks)
    Tickslabel{i} = num2str(Ticks(i));
end
cbar2.Ticks = Ticks;
cbar2.Label.String = '$log_{10}\left(\Lambda_0\right)$';
cbar2.Label.Interpreter = 'latex';
cbar2.TickLabelInterpreter = 'latex';
cbar2.Units = 'centimeters';
cbar2.Color    = [0,0,0];
%cbar2.TicksLabels = Tickslabel; 
figure_name = 'Figure_2.png';
filename = fullfile(path_save,figure_name);

exportgraphics(f,filename,'Resolution',1000,'BackgroundColor','white')         









function [ax] = set_default_Lambda(ax)
ax.XScale = 'log';
ax.XColor = [0,0,0];
ax.XLim = [10^(-6),50];
ax.XTick = [10^(-4),10^0];
ax.XTickLabel = {'$10^{-4}$','$10^{0}$'};
ax.XLabel.String = '$\Lambda_0$';
ax.XLabel.Position(2)=-0.10;
ax.XLabel.Interpreter = 'latex';
ax.LineWidth = 1.2;
ax.YLim = [0,1.1];
ax.YColor = [0,0,0]; 
ax.TickLabelInterpreter = 'latex';
ax.YLabel.String = '$1/t^{\dagger}_d$'; 

ax.YLabel.Interpreter = 'latex';
ax.YTick = [0.1,0.5,1.0]; 
ax.YTickLabel = {'$0.1$',[],'$1.0$'};
ax.Box ='on';
ax.FontUnits = 'centimeters';
ax.FontSize = 0.5;


end

function [ax] = set_default_necking(ax)
ax.XScale = 'linear';
ax.XColor = [0,0,0];
ax.XLim = [0.0,50];
ax.XTick = [0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40,45,50.0];
ax.XTickLabel = {'0',[], [],[],'$20$',[],[],[], '$40$',[],[]};
ax.XLabel.String = '$t^{\dagger}$';
ax.XTickLabelRotation = 0;
ax.XLabel.Position(2)=-0.03;
ax.XLabel.Interpreter = 'latex';
ax.LineWidth = 1.2;
ax.YLim = [0,1.1];
ax.YColor = [0,0,0]; 
ax.TickLabelInterpreter = 'latex';
ax.YLabel.String = '$D^{\dagger}$'; 
ax.YLabel.Interpreter = 'latex';
ax.YTick = [0.1,0.5,1.0]; 
ax.YLim = [0.09,1.05];
ax.YTickLabel = {'$0.1$',[],'$1.0$'};
ax.Box ='on';
ax.FontUnits = 'centimeters';
ax.FontSize = 0.5;


end

