% Figure Manuscript Non Linear portion
clear all
close all
addpath('../Utilities/')
addpath('../Utilities/Plot_Class/')
addpath('../Utilities/ScientificColourMaps8/')
clf; 
%%NonLinear_Tests_Data_Base_n3_iteration_high_T.mat It=load('NonLinear_Tests_Data_Base_n3_iteration.mat');/NonLinear_Tests_Data_Base_delta_L0
%It=load('..\Data_Base\NonLinear_Tests_Data_Base_n3_iteration_high_T_initial_guess.mat');

A = '../../Data_Base/Linear_Regime_Main_data_base_gamma_1000e3_exp.mat';
%B = '../../Data_Base/Linear_Regime_Main_data_base_gamma_1000e3.mat';
%C = '../../Data_Base/Transition_Regime_Main_data_base_gamma_1000e3.mat';
L = load(A);
%L = load(B);
%T = load(C);
S1 = L.Data_S;
%S2 = T.Data_S;
%S3 = NL.Data_S;
%%
f=figure(1)
clf;

set(gcf, 'Units','centimeters', 'Position', [0, 0, 15,15], 'PaperUnits', 'centimeters', 'PaperSize', [15, 15])
p0x = 0.08;
p0y = 0.09; 
sz_xh = 0.4;
sz_yv = 0.4; 
dh    = 0.05; 
dv    = 0.05; 

x_n = S1.Lambda0(S1.xiUM<1e-6 & S1.xiUS>1e4);
y_n = S1.tdet(S1.xiUM<1e-6 & S1.xiUS>1e4); 
N   = S1.n(S1.xiUM<1e-6 & S1.xiUS>1e4);
n   = unique(S1.n);
n_list = {[ 0.9686    0.8627    0.4353],[ 88 214 141]./255,[211 84 0]./255};
x_s = S1.Lambda0(S1.xiUM<1e-6 & S1.n>3.5);
y_s = S1.tdet(S1.xiUM<1e-6 & S1.n>3.5);
ss   = unique(S1.xiUS);
s_list = {[171 235 198]./255,[ 214 234 248 ]./255,[205 97 85]./255,[ 203 67 53]./255};
SS   = S1.xiUS(S1.xiUM<1e-6 & S1.n==3.5);




subplot(2,2,1)
% Exponents Column
ax1 = gca;
ax1.Position = [p0x,p0y+sz_yv+dv,sz_xh,sz_yv];
hold on 
for i=1:length(n)
    p(i) = scatter(x_n(N==n(i)),y_n(N==n(i)),40,n_list{i},'filled','o','MarkerEdgeColor','k');
end
hold off
ax1.XScale = 'log'; 
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1.XLabel.String = [];
ax1.XTick = [10^-5,10^-3,10^0];
ax1.XTickLabel = [];
ax1.TickLabelInterpreter = 'latex'; 
ax1.YLim    = [0,20];
ax1.YLabel.String = '$t^{\dagger}_{det}$';
ax1.YLabel.Interpreter = 'latex'; 
ax1.Box = 'on';
ax1.LineWidth =1.2;
ax1.XLim = [10^(-6),10^(0.5)]; 
l = legend('$n=2.0$','$n=3.5$','$n=7.0$');
l.Interpreter = 'latex';
l.Box = 'on';
l.Location ='west'; 
l.Position(2) = 0.85; 


% Slab dominance Column
subplot(2,2,2)
ax2 = gca;
ax2.Position = [p0x+sz_xh+dh,p0y+sz_yv+dv,sz_xh,sz_yv];
hold on 
for i=1:length(ss)
    p2(i) = scatter(x_s(SS==ss(i)),y_s(SS==ss(i)),40,s_list{i},'filled','o','MarkerEdgeColor','k');
end
hold off
ax2.XScale = 'log'; 
ax2.XColor = 'k';
ax2.YColor = 'k';
ax2.XLabel.String = [];
ax2.XTick = [10^-5,10^-3,10^0];
ax2.XTickLabel = [];
ax2.TickLabelInterpreter = 'latex'; 
ax2.YLabel.String = [];
ax2.YTickLabel = []; 
ax2.YLim    = [0,20];
ax2.Box = 'on';
ax2.LineWidth =1.2;
ax2.XLim = [10^(-6),10^(0.5)]; 
l2 = legend('$log_{10}\left(\xi_S\right)=0.0$','$log_{10}\left(\xi_S\right)=1.0$','$log_{10}\left(\xi_S\right)=2.0$','$log_{10}\left(\xi_S\right)=6.0$');
l2.Interpreter = 'latex';
l2.Box = 'on';
l2.Location ='west'; 
l2.Position(2) = 0.82


% Exponent Column
subplot(2,2,3)
ax3 = gca;
ax3.Position = [p0x,p0y,sz_xh,sz_yv];
hold on 
for i=1:length(n)
    p3(i) = scatter(x_n(N==n(i)),1./y_n(N==n(i)),40,n_list{i},'filled','o','MarkerEdgeColor','k');
end
hold off
ax3.XScale = 'log'; 
ax3.XColor = 'k';
ax3.YColor = 'k';
ax3.YLabel.String = '$\frac{1}{t^{\dagger}_{det}}$';
ax3.YLabel.Interpreter = 'latex'; 
ax3.XLabel.String = '$\Lambda_0$';
ax3.XLabel.Interpreter = 'latex'; 
ax3.XTick = [10^-5,10^-3,10^0];
ax3.XTickLabel = {'$10^{-5}$','$10^{-3}$','$10^{0}$'};
ax3.TickLabelInterpreter = 'latex'; 
ax3.YLim    = [0,1.1];
ax3.XLim = [10^(-6),10^(0.5)]; 

ax3.Box = 'on';
ax3.LineWidth =1.2;

% Exponent column
subplot(2,2,4)
ax4 = gca;
ax4.Position = [p0x+sz_xh+dh,p0y,sz_xh,sz_yv];

hold on 
for i=1:length(ss)
    p4(i) = scatter(x_s(SS==ss(i)),1./y_s(SS==ss(i)),40,s_list{i},'filled','o','MarkerEdgeColor','k');
end
hold off
ax4.XScale = 'log'; 
ax4.XColor = 'k';
ax4.YColor = 'k';
ax4.XLabel.String = '$\Lambda_0$';
ax4.XLabel.Interpreter = 'latex';
ax4.XLim = [10^(-6),10^(0.5)]; 
ax4.XTick = [10^-5,10^-3,10^0];
ax4.XTickLabel = {'$10^{-5}$','$10^{-3}$','$10^{0}$'};
ax4.TickLabelInterpreter = 'latex'; 
ax4.YLabel.String = [];
ax4.YTickLabel = []; 
ax4.YLim    = [0,1.1];
ax4.Box = 'on';
ax4.LineWidth =1.2;





t0 = annotation("textbox");%ax1.Position = [p0x,p0y+sz_yv+dv,sz_xh,sz_yv];
t0.Position = [p0x+sz_xh-0.05,p0y+2*sz_yv-0.01,0.05,0.05];
t0.Interpreter = 'latex';
t0.String = '$[a]$';
t0.FontSize = 12;
t0.FontWeight = 'bold';

t0.LineStyle = 'none';

t1 = annotation("textbox");
t1.Position = [p0x+2*sz_xh+dv-0.05,p0y+2*sz_yv-0.01,0.05,0.05];
t1.Interpreter = 'latex';
t1.String = '$[b]$';
t1.FontSize = 12;
t1.FontWeight = 'bold';
t1.LineStyle = 'none';

t2 = annotation("textbox");%ax1.Position = [p0x,p0y+sz_yv+dv,sz_xh,sz_yv];
t2.Position = [p0x+sz_xh-0.05,p0y+1*sz_yv-0.05,0.05,0.05];
t2.Interpreter = 'latex';
t2.String = '$[c]$';
t2.FontSize = 12;
t2.FontWeight = 'bold';

t2.LineStyle = 'none';

t3 = annotation("textbox");
t3.Position = [p0x+2*sz_xh+dv-0.05,p0y+1*sz_yv-0.05,0.05,0.05];
t3.Interpreter = 'latex';
t3.String = '$[d]$';
t3.FontSize = 12;
t3.FontWeight = 'bold';
t3.LineStyle = 'none';



exportgraphics(f,'Figure_S2.png','Resolution',1000,'BackgroundColor','white')         








function [ax] = set_default_Lambda(ax)
ax.XScale = 'log';
ax.XColor = [0,0,0];
ax.XLim = [10^(-6),50];
ax.XTick = [10^(-4),10^0];
ax.XTickLabel = {'$10^{-4}$','$10^{0}$'};
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

