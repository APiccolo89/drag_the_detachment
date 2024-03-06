% Figure Manuscript Non Linear portion
clear all
close all
addpath('../Utilities/')
addpath('../Utilities/Plot_Class/')
addpath('../Utilities/ScientificColourMaps8/')
clf; 
%%NonLinear_Tests_Data_Base_n3_iteration_high_T.mat It=load('NonLinear_Tests_Data_Base_n3_iteration.mat');/NonLinear_Tests_Data_Base_delta_L0
%It=load('..\Data_Base\NonLinear_Tests_Data_Base_n3_iteration_high_T_initial_guess.mat');

C = '../../Data_Base/Transition_Regime_Main_data_base_gamma_1000e3.mat';
T = load(C);
S2 = T.Data_S;
fun_0D = Manuscript_function_Container;  

[Transition_tau_M]=fun_0D.select_tests_prepare_variables(T.n_3,0,'time_nd','tau_M','Lambda0','NonLinear',3,'xius',[]);
[Transition_partition]=fun_0D.select_tests_prepare_variables(T.n_3,0,'time_nd','xdisl','Lambda0','NonLinear',3,'xius',[]);


%%  Lambda Figure 
%%  Lambda Figure


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
Tcolor = fun_0D.color_computation(3,transitional_critical_lambda,-4,1);
path2colormap = strcat('..\Utilities\ScientificColourMaps8\','lipari','\','lipari','.mat');

load(path2colormap);
cmap        = colormap(lipari);

close all;

f=figure(1);
clf;
set(gcf, 'Units','centimeters', 'Position', [0, 0, 15,15], 'PaperUnits', 'centimeters', 'PaperSize', [15, 15])
% linear experiment
subplot(1,2,1)
% Function that creates the data structure with the relevant information
% for the picture at hand 
path2colormap = strcat(['..\Utilities\ScientificColourMaps8\'],'lipari','\','lipari','.mat');

load(path2colormap);

cmap        = colormap(lipari);

size_tests = length(squeeze(Transition_tau_M(1,1,:)));
c = squeeze(Transition_tau_M(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c(1,:),-4,1);
subplot(2,3,4)
ax3 = gca(); 
ax3.Position = [0.08,0.1,0.40,0.8];
x = squeeze(Transition_tau_M(1,:,:));
y = squeeze(Transition_tau_M(2,:,:));
c = squeeze(Transition_tau_M(3,:,:));
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
       if abs(log10(cc(1))-log10(transitional_critical_lambda(1)))<1e-2 || abs(log10(cc(1))-log10(transitional_critical_lambda(2)))<1e-3 || abs(log10(cc(1))-log10(transitional_critical_lambda(3)))<1e-3
            disp(cc(1));
            pl2(i)=plot(ax3,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',1.2);
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





subplot(1,2,2)
ax4 = gca(); 
ax4.Position = [0.56,0.1,0.40,0.8];


size_tests = length(squeeze(Transition_partition(1,1,:)));
c = squeeze(Transition_partition(3,:,:));
color_lists = fun_0D.color_computation(size_tests,c(1,:),-4,1);

x = squeeze(Transition_partition(1,:,:));
y = squeeze(Transition_partition(2,:,:));
c = squeeze(Transition_partition(3,:,:));



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

ax4 = set_default_partition(ax4); 

ax4.YTickLabel = [];
ax4.YLabel.String = [];
caxis(ax4,[-4,1])
yline(ax4,1.0,'LineWidth',1.2,'color','k','Alpha',0.3);
xline(ax4,1.0,'LineWidth',1.2,'Color','b','Alpha',0.3)




cbar2= colorbar(ax4,"southoutside");
cbar2.Limits = [-4 1];
cbar2.Position = [0.25,0.94,0.50,0.015];

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


print(f,filename,'-dpng','-r0')






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
ax.YLabel.String = '$\tau^{\dagger}_M$'; 
ax.YLabel.Interpreter = 'latex';
ax.YScale = 'log';
ax.YTick = [0.0,0.5,1.0]; 
ax.YLim = [0.0,1.0];
ax.YTickLabel = {'$0.0$','0.5','$1.0$'};
ax.Box ='on';
ax.FontUnits = 'centimeters';
ax.FontSize = 0.5;

end

function [ax] = set_default_partition(ax)
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
ax.YLabel.String = '$\nu_{disl}$'; 
ax.YLabel.Interpreter = 'latex';
ax.YScale = 'log';
ax.YTick = [0.0,0.5,1.0]; 
ax.YLim = [0.001,1.0];
ax.YTickLabel = {'$0.0$','$5.0$','$10.0$'};
ax.Box ='on';
ax.FontUnits = 'centimeters';
ax.FontSize = 0.5;


end

