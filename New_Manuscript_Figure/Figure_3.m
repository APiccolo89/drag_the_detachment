% Figure Manuscript Non Linear portion
clear all
close all
addpath('../Utilities/')
addpath('../Utilities/Plot_Class/')
addpath('../Utilities/ScientificColourMaps8/')
clf; 
%%NonLinear_Tests_Data_Base_n3_iteration_high_T.mat It=load('NonLinear_Tests_Data_Base_n3_iteration.mat');/NonLinear_Tests_Data_Base_delta_L0
%It=load('..\Data_Base\NonLinear_Tests_Data_Base_n3_iteration_high_T_initial_guess.mat');

A = '../../Data_Base/NonLinear_Regime_Main_data_base_gamma_1000e3.mat';
B = '../../Data_Base/Linear_Regime_Main_data_base_gamma_1000e3.mat';
C = '../../Data_Base/Transition_Regime_Main_data_base_gamma_1000e3.mat';
NL = load(A);
L = load(B);
T = load(C);
S1 = L.Data_S;
S2 = T.Data_S;
S3 = NL.Data_S;
fun_0D = Manuscript_function_Container;  

[Linear]=fun_0D.select_tests_prepare_variables(L.n_3,0,'time_nd','tau_eff','Lambda0','NonLinear',3,'xius',[]);
[Transition]=fun_0D.select_tests_prepare_variables(T.n_3,0,'time_nd','tau_eff','Lambda0','NonLinear',3,'xius',[]);
[Nonlinear]=fun_0D.select_tests_prepare_variables(NL.n_3,0,'time_nd','tau_eff','Lambda0','NonLinear',3,'xius',[]);
%%  Lambda Figure 
%%  Lambda Figure
% Select a few data per each group:
buf1 = find(1./S1.tdet>0.95);
buf2 = median(log10(S1.Lambda0(buf1)));
buf2 = round(buf2);
buf3 =find(1./S1.tdet>0.45 & 1./S1.tdet<0.55);
buf4 = median(log10(S1.Lambda0(buf3)));
buf4 = round(buf4,3);
buf5 =find(1./S1.tdet>0.05 & 1./S1.tdet<0.15);
buf6 = median(log10(S1.Lambda0(buf5)));
buf6 = round(buf6,3);
linear_critical_lambda = 10.^[buf2,buf4,buf6];


buf1 = find(1./S3.tdet>0.95);
buf2 = median(log10(S3.Lambda0(buf1)));
buf2 = round(buf2,3);

buf3 =find(1./S3.tdet>0.45 & 1./S3.tdet<0.55);
buf4 = median(log10(S3.Lambda0(buf3)));
buf4 = round(buf4,3);

buf5 =find(1./S3.tdet>0.05 & 1./S3.tdet<0.15);
buf6 = median(log10(S3.Lambda0(buf5)));
buf6 = round(buf6,3);

Nlinear_critical_lambda = 10.^[buf2,buf4,buf6];


buf1 = find(1./S2.tdet>0.95);
buf2 = median(log10(S2.Lambda0(buf1)));
buf2 = round(buf2,3);

buf3 =find(1./S2.tdet>0.45 & 1./S2.tdet<0.55);
buf4 = median(log10(S2.Lambda0(buf3)));
buf4 = round(buf4,3);

buf5 =find(1./S2.tdet>0.05 & 1./S2.tdet<0.15);
buf6 = median(log10(S2.Lambda0(buf5)));
buf6 = round(buf6,3);


transitional_critical_lambda = 10.^[buf2,buf4,buf6];
Lcolor = fun_0D.color_computation(3,linear_critical_lambda,-4,1);
Tcolor = fun_0D.color_computation(3,transitional_critical_lambda,-4,1);
NLcolor = fun_0D.color_computation(3,Nlinear_critical_lambda,-4,1);
path2colormap = strcat('..\Utilities\ScientificColourMaps8\','lipari','\','lipari','.mat');

load(path2colormap);
cmap        = colormap(lipari);

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

s1 = scatter(ax0,S1.Lambda0,S1.tau_max,30,log10(S1.xiUM),"filled","MarkerEdgeColor",'k'); 
yline(1.0,'LineWidth',1.2,'Color','r','LineStyle',':')
xline(linear_critical_lambda(1),'LineWidth',1.1,'Color',cmap(Lcolor(1),:),'LineStyle',':')
xline(linear_critical_lambda(2),'LineWidth',1.1,'Color',cmap(Lcolor(2),:),'LineStyle',':')
xline(linear_critical_lambda(3),'LineWidth',1.1,'Color',cmap(Lcolor(3),:),'LineStyle',':')

[ax0] = set_default_Lambda(ax0);
ax0.YLabel.Units = 'normalized';
ax0.YLabel.Position(1)=-0.05;
ax0.YLabel.Position(2)=0.5;% This field is automatically set as a function of the DATA 
caxis(ax0,[-12,12]);
colormap(ax0,crameri('oslo'));
%l = line(ax0,[1e-8,10^max_l0],[0.5,0.50],'Color','k','LineWidth',1.2);
%str = {['I'],['n'],['v'],['i'],['s'],['c'],['i'],['d']};
t_type = text(0.9,0.1,'Inviscid','Interpreter','latex','rotation',90);
%t_type.Position(1:2) = ax0.Position(1:2); 
%t_type.Position(2) = 0.68;
%t_type.Position(1) = ax0.Position(1)+0.01;
t_type.LineStyle = 'none';
t_type.Units = 'normalized';
t_type.Position(1)=0.1;
t_type.Position(2)=0.40;
%t_type.BackgroundColor = 'w';


%l2 = line(ax0,[10^max_l0,50],[1.0,1.0],'Color','k','LineWidth',1.2);
t_type2 = text(0.9,0.1,'Drag dominated','Interpreter','latex','rotation',90);
%t_type.Position(1:2) = ax0.Position(1:2); 
%t_type.Position(2) = 0.68;
%t_type.Position(1) = ax0.Position(1)+0.01;
t_type2.LineStyle = 'none';
t_type2.Units = 'normalized';
t_type2.Position(1)=0.9;
t_type2.Position(2)=0.25;

t_type3 = text(0.25,0.1,'Linear, $\xi_M \leq 10^{-4}$','Interpreter','latex');
t_type3.LineStyle = 'none';
t_type3.Units = 'normalized';
t_type3.Position(1)=0.15;
t_type3.Position(2)=1.1;
subplot(2,3,2)
% lNinear experiment
ax1 = gca();
p2 = ax0.Position(1)+ax0.Position(3)+0.025;
ax1.Position = [p2,0.65,0.28,0.3];

s2 = scatter(ax1,S2.Lambda0,S2.tau_max,30,log10(S2.xiUM),"filled","MarkerEdgeColor",'k'); 
yline(1.0,'LineWidth',1.2,'Color','r','LineStyle',':')
xline(transitional_critical_lambda(1),'LineWidth',1.1,'Color',cmap(Tcolor(1),:),'LineStyle',':')
xline(transitional_critical_lambda(2),'LineWidth',1.1,'Color',cmap(Tcolor(2),:),'LineStyle',':')
xline(transitional_critical_lambda(3),'LineWidth',1.1,'Color',cmap(Tcolor(3),:),'LineStyle',':')
[ax1] = set_default_Lambda(ax1);
ax1.YTickLabel = [];
ax1.YLabel.String = [];
caxis(ax1,[-12,12]);
colormap(ax1,crameri('oslo'));
t_type01 = text(0.9,0.1,'Inviscid','Interpreter','latex','rotation',90);
%t_type.Position(1:2) = ax0.Position(1:2); 
%t_type.Position(2) = 0.68;
%t_type.Position(1) = ax0.Position(1)+0.01;
t_type01.LineStyle = 'none';
t_type01.Units = 'normalized';
t_type01.Position(1)=0.1;
t_type01.Position(2)=0.40;
%t_type.BackgroundColor = 'w';


%l2 = line(ax0,[10^max_l0,50],[1.0,1.0],'Color','k','LineWidth',1.2);
t_type02 = text(0.9,0.1,'Drag dominated','Interpreter','latex','rotation',90);
%t_type.Position(1:2) = ax0.Position(1:2); 
%t_type.Position(2) = 0.68;
%t_type.Position(1) = ax0.Position(1)+0.01;
t_type02.LineStyle = 'none';
t_type02.Units = 'normalized';
t_type02.Position(1)=0.9;
t_type02.Position(2)=0.25;

t_type03 = text(0.25,0.1,'Transitional, $10^{-4}<\xi_M < 10^{4}$','Interpreter','latex');
t_type03.LineStyle = 'none';
t_type03.Units = 'normalized';
t_type03.Position(1)=-0.03;
t_type03.Position(2)=1.1;
subplot(2,3,3)
% lNinear experiment
ax2 = gca();
p3 = ax1.Position(1)+ax1.Position(3)+0.025;

ax2.Position = [p3,0.65,0.28,0.3];
s3 = scatter(ax2,S3.Lambda0,S3.tau_max,30,log10(S3.xiUM),"filled","MarkerEdgeColor",'k'); 
yline(1.0,'LineWidth',1.2,'Color','r','LineStyle',':')
xline(linear_critical_lambda(1),'LineWidth',1.1,'Color',cmap(NLcolor(1),:),'LineStyle',':')
xline(linear_critical_lambda(2),'LineWidth',1.1,'Color',cmap(NLcolor(2),:),'LineStyle',':')
xline(linear_critical_lambda(3),'LineWidth',1.1,'Color',cmap(NLcolor(3),:),'LineStyle',':')
[ax2] = set_default_Lambda(ax2);
ax2.YTickLabel = [];
ax2.YLabel.String = [];
caxis(ax2,[-12,12]);
colormap(ax2,crameri('oslo'));

t_type11 = text(0.9,0.1,'Inviscid','Interpreter','latex','rotation',90);
%t_type.Position(1:2) = ax0.Position(1:2); 
%t_type.Position(2) = 0.68;
%t_type.Position(1) = ax0.Position(1)+0.01;
t_type11.LineStyle = 'none';
t_type11.Units = 'normalized';
t_type11.Position(1)=0.1;
t_type11.Position(2)=0.40;
%t_type.BackgroundColor = 'w';


%l2 = line(ax0,[10^max_l0,50],[1.0,1.0],'Color','k','LineWidth',1.2);
t_type12 = text(0.9,0.1,'Drag dominated','Interpreter','latex','rotation',90);
%t_type.Position(1:2) = ax0.Position(1:2); 
%t_type.Position(2) = 0.68;
%t_type.Position(1) = ax0.Position(1)+0.01;
t_type12.LineStyle = 'none';
t_type12.Units = 'normalized';
t_type12.Position(1)=0.9;
t_type12.Position(2)=0.25;

t_type13 = text(0.25,0.1,'Nonlinear, $\xi_M \geq 10^{4}$','Interpreter','latex');
t_type13.LineStyle = 'none';
t_type13.Units = 'normalized';
t_type13.Position(1)=0.15;
t_type13.Position(2)=1.1;
% Function that creates the data structure with the relevant information
% for the picture at hand 
path2colormap = strcat(['..\Utilities\ScientificColourMaps8\'],'lipari','\','lipari','.mat');

load(path2colormap);

cmap        = colormap(lipari);

size_tests = length(squeeze(Linear(1,1,:)));
c = squeeze(Linear(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c(1,:),-4,1);
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
        pl(i)=plot(ax3,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.4,LineStyle=':');
        if abs(log10(cc(1))-log10(linear_critical_lambda(1)))<1e-2 || abs(log10(cc(1))-log10(linear_critical_lambda(2)))<1e-3 || abs(log10(cc(1))-log10(linear_critical_lambda(3)))<1e-3
            disp(cc(1));
            if log10(cc(1))==log10(linear_critical_lambda(2))
                bla = 0; 
            end
            pl(i)=plot(ax3,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',1.2);
        end
    end
    hold off
end
colormap(ax3,cmap);
ax3 = set_default_necking(ax3); 
ax3.YLabel.Units = 'normalized';
ax3.YLabel.Position(1)=0.0000005;
ax3.YLabel.Position(2)=0.5;
caxis(ax3,[-4,1])
yline(ax3,1.0,'LineWidth',1.2,'color','k','Alpha',0.3);
xline(ax3,1.0,'LineWidth',1.2,'Color','b','Alpha',0.3)





subplot(2,3,5)
ax4 = gca(); 
ax4.Position = [p2,0.2,0.28,0.35];


size_tests = length(squeeze(Transition(1,1,:)));
c = squeeze(Transition(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c(1,:),-4,1);

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
        pl22(i)=plot(ax4,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.4,LineStyle=':');
         if abs(log10(cc(1))-log10(transitional_critical_lambda(1)))<1e-2 || abs(log10(cc(1))-log10(transitional_critical_lambda(2)))<1e-3 || abs(log10(cc(1))-log10(transitional_critical_lambda(3)))<1e-3
            disp(cc(1));
            pl22(i)=plot(ax4,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',1.2);
        end
    end
    hold off
end

ax4 = set_default_necking(ax4); 

ax4.YTickLabel = [];
ax4.YLabel.String = [];
caxis(ax4,[-4,1])
yline(ax4,1.0,'LineWidth',1.2,'color','k','Alpha',0.3);
xline(ax4,1.0,'LineWidth',1.2,'Color','b','Alpha',0.3)


subplot(2,3,6)
ax5 = gca(); 
ax5.Position = [p3,0.2,0.28,0.35];


size_tests = length(squeeze(Nonlinear(1,1,:)));
c = squeeze(Nonlinear(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c(1,:),-4,1);

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
        pl33(i)=plot(ax5,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.4,LineStyle=':');
        if abs(log10(cc(1))-log10(Nlinear_critical_lambda(1)))<1e-2 || abs(log10(cc(1))-log10(Nlinear_critical_lambda(2)))<1e-3 || abs(log10(cc(1))-log10(Nlinear_critical_lambda(3)))<1e-3            disp(cc(1));
            pl33(i)=plot(ax5,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',1.2);
        end
    end
    hold off
end

ax5 = set_default_necking(ax5); 
ax5.YTickLabel = [];
ax5.YLabel.String = [];
caxis(ax5,[-4,1])
yline(ax5,1.0,'LineWidth',1.2,'color','k','Alpha',0.3);
xline(ax5,1.0,'LineWidth',1.2,'Color','b','Alpha',0.3)


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
figure_name = 'Figure_3.png';
filename = fullfile(figure_name);



t0 = annotation("textbox");
t0.Position = [0.08+0.22,0.65+0.24,0.05,0.05];
t0.Interpreter = 'latex';
t0.String = '$[a]$';
t0.FontSize = 12;
t0.FontWeight = 'bold';

t0.LineStyle = 'none';

t1 = annotation("textbox");
t1.Position = [0.385+0.22,0.65+0.24,0.05,0.05];
t1.Interpreter = 'latex';
t1.String = '$[b]$';
t1.FontSize = 12;
t1.FontWeight = 'bold';
t1.LineStyle = 'none';


t2= annotation("textbox");
t2.Position = [0.69+0.22,0.65+0.24,0.05,0.05];
t2.Interpreter = 'latex';
t2.String = '$[c]$';
t2.FontSize = 12;
t2.FontWeight = 'bold';
t2.LineStyle = 'none';



t3 = annotation("textbox");
t3.Position = [0.08+0.22,0.65-0.15,0.05,0.05];%0.08+0.22
t3.Interpreter = 'latex';
t3.String = '$[d]$';
t3.FontSize = 12;
t3.FontWeight = 'bold';

t3.LineStyle = 'none';

t4 = annotation("textbox");
t4.Position = [0.385+0.22,0.65-0.15,0.05,0.05];
t4.Interpreter = 'latex';
t4.String = '$[e]$';
t4.FontSize = 12;
t4.FontWeight = 'bold';
t4.LineStyle = 'none';


t5= annotation("textbox");
t5.Position = [0.69+0.22,0.65-0.15,0.05,0.05];
t5.Interpreter = 'latex';
t5.String = '$[f]$';
t5.FontSize = 12;
t5.FontWeight = 'bold';
t5.LineStyle = 'none';

exportgraphics(f,filename,'Resolution',1000,'BackgroundColor','white')         








function [ax] = set_default_Lambda(ax)
ax.XScale = 'log';
ax.XColor = [0,0,0];
ax.XLim = [10^(-8),50];
ax.XTick = [10^(-7),10^(-4),10^0];
ax.XTickLabel = {'$10^{-7}$','$10^{-4}$','$10^{0}$'};
ax.XLabel.String = '$\Lambda_0$';

ax.XLabel.Units = 'normalized';
ax.XLabel.Position(2) = -0.10;

ax.XLabel.Interpreter = 'latex';
ax.LineWidth = 1.2;
ax.YLim = [0,11];
ax.YColor = [0,0,0]; 
ax.TickLabelInterpreter = 'latex';
ax.YLabel.String = '$\tau^{\dagger}_{MAX}$'; 

ax.YLabel.Interpreter = 'latex';
ax.YTick = [0.0, 0.5,1.0,5.0,10.0]; 
ax.YTickLabel = {'$0.0$',[],'$1.0$',[],'$10.0$'};
ax.Box ='on';
ax.FontUnits = 'centimeters';
ax.FontSize = 0.5;

end

function [ax] = set_default_necking(ax)
ax.XScale = 'log';
ax.XColor = [0,0,0];
ax.XLim = [0.1,50];
ax.XTick = [0.1000    1.0000   10.0000   40.0000];
ax.XTickLabel = {'0.1','1','10',[]};
ax.XLabel.String = '$t^{\dagger}$';
ax.XTickLabelRotation = 0;
ax.XLabel.Units = 'normalized';
ax.XLabel.Position(2)=-0.1;
ax.XLabel.Interpreter = 'latex';
ax.LineWidth = 1.2;
ax.YColor = [0,0,0]; 
ax.TickLabelInterpreter = 'latex';
ax.YLabel.String = '$\tau^{\dagger}$'; 
ax.YLabel.Interpreter = 'latex';
ax.YScale = 'log';
ax.YTick = [0.1,1.0,5.0,10.0]; 
ax.YLim = [0.5,11.0];
ax.YTickLabel = {'$0.0$','1.0','5.0','$10.0$'};
ax.Box ='on';
ax.FontUnits = 'centimeters';
ax.FontSize = 0.5;


end

