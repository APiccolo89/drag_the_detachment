clear all;
close all;
addpath ../Utilities/
addpath ../../../../../freezeColors-master/



% Comparison 0D - 2D
load('../Data_Base_Fit.mat')

path_save = '../../Manuscript2';


%plot_comparison(NLN.TB.T1,NLN2.TB.T1,path_save,'Comparison_T1',[0,3.0])
%plot_comparison(NLN.TB.T15,NLN2.TB.T15,path_save,'Comparison_T15',[0,9.0])
% Plot data Linear Experiment database
% create path where to save if does not exit
nlm = Problem_type.Linear;
pt = fullfile(path_save,'Maps_NLinear');
if ~isfolder(pt)
    mkdir(pt)
end
%%
% Create figure
%plot_the_figures(LN_ALT,pt,nlm)
%[linear]=plot_2Dmaps_profiles(NLN_ALT,pt,nlm);
%[Nlinear]=plot_2Dmaps_profiles(LN_ALT,pt,nlm);
% %%
% figure(2)
% 
% clf;
% ax = gca; 
% Data = [ones(length(linear),1)'];
% Data2 = [ones(length(Nlinear),1)'.*2];
% hold on; 
% boxchart(Data,linear,"BoxFaceColor",[ 230, 126, 34]./255,'Notch','on')
% boxchart(Data2,Nlinear,"BoxFaceColor",[ 33, 97, 140]./255,'Notch','on')
% 
% ax.Box = 'on';
% ax.LineWidth = 1.2; 
% ax.XGrid = 'on';
% ax.YGrid = 'on';
% ax.XColor = [0 0 0 ];
% ax.YColor = [0 0 0 ]; 
% figure_name = 'box_plot.png';
% path_2_save = fullfile(path_save,figure_name);
% print(path_2_save,'-dpng','-r600')
% 



% % Plot data Non Linear Experiment with one fetches
 pt = fullfile(path_save,'NonLinear_Comparison_New');
 nlm = Problem_type.NonLinear;
% 
if ~isfolder(pt)
    mkdir(pt)
end
plot_the_figures(NLN,pt,nlm)

% % Plot data Non linear experiments with two fetches
% pt = fullfile(path_save,'NonLinear_Comparison_2fetches_New');
% nlm = Problem_type.NonLinear;
% 
% 
% if ~isfolder(pt)
%     mkdir(pt)
% end
% 
% plot_the_figures(LN_ALT,pt,nlm)

pt = fullfile(path_save,'Linear_Comparison_2fetches_New');
nlm = Problem_type.Linear;


if ~isfolder(pt)
    mkdir(pt)
end

plot_the_figures(LN,pt,nlm)

% Figure Comparison Non linear experiment

function plot_comparison(Tests,Test2fetches,pt_save,name_figure,tlim)


clf;
close all;

font_axes = 16;
font_legend = 14;
font_text   = 5;
size_shit = [12,13.5];
LineWidth = 1.0;
marker_size = 10;
F1 = figure(1);

t_d = 0:0.001:1;
D   = (1-t_d).^(1/3.5);

set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])


plot(Tests.t.*3.5,Tests.D,"Color",'#9a2a3c',LineWidth=1.2,LineStyle="-")
hold on
plot(Tests.t.*3.5,Tests.D0D2,'Color','#b2ec73',LineWidth=1.0,LineStyle='--')
plot(Test2fetches.t.*3.5,Test2fetches.D0D2,'Color','#6e4e8e',LineWidth=1.0,LineStyle='--')
plot(t_d,D,'Color','k','LineStyle','-.','LineWidth',1.0)

ax2 = gca;
ax2.XLim = tlim;
ax2.YLim = [0.1,1.1];
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.FontSize = font_axes;
ax2.TickLabelInterpreter = 'latex';
ax2.Box  = 'on';
ax2.XGrid = 'on';
ax2.YGrid = 'on';
ax2.XColor = [0 0 0 ];
ax2.YColor = [0 0 0 ];
ax2.LineWidth = 1.2;
l=legend('$2D$','$0D(\omega_{depth})$','$0D(\omega_S,\omega_{depth})$','Analytical solution');
l.Interpreter = 'latex';
l.FontSize = font_legend;
l.Location = 'northeast';

pt=fullfile(pt_save,name_figure);
print(pt,'-dpng')



end


function plot_the_figures(S,pt_save,nlm)
close all;
clf;
FIT = S.FIT;
TB  = S.TB;
font_axes = 16;
font_legend = 14;
font_text   = 5;
size_shit = [12,13.5];
LineWidth = 1.2;
marker_size = 10;

% figure(1)
% clf;
% set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
% ax = gca;
% hold on
% scatter(FIT.Lambda,1-FIT.fetch(1,:),20,[204 255 204]./255,'filled','d','MarkerEdgeColor','k')
% ax.XScale = 'log';
% ax.XColor = [0,0,0];
% ax.XLim = [10^-6,10^0];
% ax.YLim = [0.1, 0.5];
% %ax.XLabel.String = "$\frac{\Lambda}{1+\xi^{UM}}$";
% ax.XLabel.Interpreter = 'latex';
% ax.YColor = [0,0,0];
% %ax.YLabel.String = "$\omega_{\mathrm{depth}}$";
% ax.YLabel.Interpreter = 'latex';
% ax.XGrid = 'on';
% %ax.XTickLabel = [];
% %ax.YTickLabel = [];
% ax.XMinorTick = 'on';
% ax.XMinorGrid = 'on';
% ax.YGrid = 'on';
% ax.YMinorTick = 'on';
% ax.YMinorGrid = 'on';
% ax.Box = 'on';
% ax.TickLabelInterpreter = 'latex';
% ax.LineWidth = LineWidth;
% ax.FontSize = font_axes;
% %l = legend('2D simulations','0D fitting');
% %l.Interpreter = 'latex';
% figure_name ='Depth_Coefficient';
% 
% pt=fullfile(pt_save,figure_name);
% print(pt,'-dpng')
% %%
% figure(2)
% clf;
% set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
% ax = gca;
% hold on
% scatter(FIT.Lambda,abs((1-FIT.fetch(1,:))./FIT.Dp./2),10,FIT.L0./1e3,'filled','d'); colorbar;
% %scatter(FIT.Lambda,(-FIT.Dp(1,:))./2,10,'k','filled','+');
% 
% ax.XScale = 'log';
% ax.XColor = [0,0,0];
% ax.XLim = [10^-6,10^0];
% %ax.YLim = [0.1, 0.5];
% %ax.XLabel.String = "$\frac{\Lambda}{1+\xi^{UM}}$";
% ax.XLabel.Interpreter = 'latex';
% ax.YColor = [0,0,0];
% %ax.YLabel.String = "$\omega_{\mathrm{depth}}$";
% ax.YLabel.Interpreter = 'latex';
% ax.XGrid = 'on';
% %ax.XTickLabel = [];
% %ax.YTickLabel = [];
% ax.XMinorTick = 'on';
% ax.XMinorGrid = 'on';
% ax.YGrid = 'on';
% ax.YMinorTick = 'on';
% ax.YMinorGrid = 'on';
% ax.Box = 'on';
% ax.TickLabelInterpreter = 'latex';
% ax.LineWidth = LineWidth;
% ax.FontSize = font_axes;
% %l = legend('2D simulations','0D fitting');
% %l.Interpreter = 'latex';
% figure_name ='Depth_Coefficient_difference';
% 
% pt=fullfile(pt_save,figure_name);
% print(pt,'-dpng')
% %%
% figure(2)
% clf;
% set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
% ax = gca;
% scatter(FIT.Lambda,FIT.fetch(2,:),20,log10(FIT.xium),'filled','o','MarkerEdgeColor','k');
% colormap(crameri('oslo'));
% hold on
% ax.XScale = 'log';
% ax.YScale = 'log';
% ax.XLim = [10^-6,10^0];
% %ax.XTickLabel = [];
% %ax.YTickLabel = [];
% ax.XColor = [0,0,0];
% %ax.YLim = [0,2];
% %ax.XLabel.String = "$\frac{\Lambda}{1+\xi^{UM}}$";
% ax.XLabel.Interpreter = 'latex';
% ax.YColor = [0,0,0];
% %ax.YLabel.String = "$\omega_{\mathrm{S}}$";
% ax.YLabel.Interpreter = 'latex';
% ax.XGrid = 'on';
% ax.XMinorTick = 'on';
% ax.XMinorGrid = 'on';
% ax.YGrid = 'on';
% ax.YMinorTick = 'on';
% ax.YMinorGrid = 'on';
% ax.YScale = 'log';
% ax.Box = 'on';
% ax.TickLabelInterpreter = 'latex';
% ax.LineWidth = LineWidth;
% ax.FontSize = font_axes;
% figure_name ='S_Coefficient';
% colorbar;
% 
% pt=fullfile(pt_save,figure_name);
% print(pt,'-dpng')
% 
% 
% 
% 
% 
% 
% figure(3)
% set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
% ax = gca;
% scatter(FIT.Lambda,FIT.fitting_p,20,[204 204 255]./255,'filled','o','MarkerEdgeColor','k');
% ax.XScale = 'log';
% ax.XColor = [0,0,0];
% %ax.XTickLabel = [];
% %ax.YTickLabel = [];
% ax.XLim = [10^-6,10^0];
% ax.YLim = [0.0, 0.4];
% ax.XLabel.String = "$\frac{\Lambda}{1+\xi^{UM}}$";
% ax.XLabel.Interpreter = 'latex';
% ax.YColor = [0,0,0];
% %ax.YLabel.String = "$F$";
% %ax.YLabel.Interpreter = 'latex';
% ax.XGrid = 'on';
% ax.XMinorTick = 'on';
% ax.XMinorGrid = 'on';
% ax.YGrid = 'on';
% ax.YMinorTick = 'on';
% ax.YMinorGrid = 'on';
% ax.Box = 'on';
% ax.TickLabelInterpreter = 'latex';
% ax.LineWidth = 1.2;
% ax.FontSize = 16;
% %l = legend('2D simulations','0D fitting');
% %l.Interpreter = 'latex';
% figure_name ='Fitting';
% 
% pt=fullfile(pt_save,figure_name);
% print(pt,'-dpng')
% 

%% 1D line plot
fnames = fieldnames(TB);
flength = length(fnames);
load('..\Utilities\ScientificColourMaps8\lipari\lipari.mat');

for ktest = 1:flength
    t = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).ID.n;
    t_dim = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).tc;
    t_dim = t_dim./(365.25*60*60*24*1e6);
    t0D = TB.(fnames{ktest}).t0D.*TB.(fnames{ktest}).ID.n;
    D = TB.(fnames{ktest}).D;
    D_0D = TB.(fnames{ktest}).D0D2;
    tau = TB.(fnames{ktest}).tau;
    tau0D = TB.(fnames{ktest}).tau0D;
    topo = TB.(fnames{ktest}).topo;
    tau_L = TB.(fnames{ktest}).tauL;
    dt_dim = t_dim(2:1:end)-t_dim(1:1:end-1);
    dTopo = topo(2:1:length(t_dim))-topo(1:1:length(t_dim)-1);
    dU = dTopo./dt_dim;
    T_Mean = (t_dim(2:1:end)+t_dim(1:end-1)).*0.5;
    t_M = 0.5.*(t(1:1:end-1)+t(2:1:end));

    if ~isempty(D_0D)
        c_l1 = [-5,0];% lim color Lambda
        c_l2 = [0,30];% lim color Res
        R    = TB.(fnames{ktest}).res.*100;
        L    = TB.(fnames{ktest}).ID.Lambda;
        if nlm.islinear == 0
            L = L./(1+TB.(fnames{ktest}).ID.Df_UM);
        end
        index_color1 = find_color(L,c_l1,1);
        index_color2 = find_color(R,c_l2,0.5);

        cmap = colormap(lipari);

%         F7 = figure(9);
%         set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%         hold on
%         colormap(crameri('bilbao',256));freezeColors;
% 
%         plot(t,tau_L(1:length(t)),"Color",cmap(index_color1,:),LineWidth=0.8,LineStyle="-")
%         ax = gca;
% 
%         ax.XLim = [0,10];
%         ax.YLim = [0,1.2];
%         % ax.XLabel.String = '$t^{\dagger}$';
%         %ax.YLabel.String = '$\frac{\tau_{lith}}{\tau_{0,B}}$';
%         ax.XLabel.Interpreter = 'latex';
%         ax.YLabel.Interpreter = 'latex';
%         ax.FontSize = font_axes;
%         ax.TickLabelInterpreter = 'latex';
%         ax.Box  = 'on';
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.XColor = [0 0 0 ];
%         ax.YColor = [0 0 0 ];
%         ax.LineWidth = 1.2;
%         % ax.XTickLabel =[];
%         % ax.YTickLabel = [];
%         ax.YScale = 'log';
%         cmap = colormap(lipari);

        % Topography
        F1 = figure(4);
        set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

        hold on
        plot(t_dim,topo(1:length(t))./(abs(min(topo(1:length(t))))),"Color",cmap(index_color1,:),LineWidth=0.8,LineStyle="-")
        ax = gca;

        %ax.XLim = [0,10];
        ax.YLim = [-1.0,0.5];
        %ax.XLabel.String = '$t^{\dagger}$';
        %ax.YLabel.String = '${H} [\mathrm{km}]$';
        ax.XLabel.Interpreter = 'latex';
        ax.YLabel.Interpreter = 'latex';
        ax.FontSize = font_axes;
        ax.TickLabelInterpreter = 'latex';
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.Box  = 'on';
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.XColor = [0 0 0 ];
        ax.YColor = [0 0 0 ];
        ax.LineWidth = 1.2;
        % ax.XTickLabel =[];
        %ax.YTickLabel = [];

%         F2 = figure(5);
%         set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%         hold on
%         plot(t_dim,topo(1:length(t)),"Color",cmap(index_color1,:),LineWidth=0.8,LineStyle="-")
%         ax2 = gca;
% 
%         %ax.XLim = [0,10];
%         %ax.YLim = [-1.5,1.5];
%         %ax2.XLabel.String = '$t [\mathrm{Myrs}]$';
%         %ax2.YLabel.String = '${H} [\mathrm{km}]$';
%         ax2.XLabel.Interpreter = 'latex';
%         ax2.YLabel.Interpreter = 'latex';
%         ax2.FontSize = font_axes;
%         ax2.TickLabelInterpreter = 'latex';
%         ax2.XMinorTick = 'on';
%         ax2.YMinorTick = 'on';
% 
%         ax2.Box  = 'on';
%         ax2.XGrid = 'on';
%         ax2.YGrid = 'on';
%         ax2.XColor = [0 0 0 ];
%         ax2.YColor = [0 0 0 ];
%         ax2.LineWidth = 1.2;
% 
%         cmap = colormap(crameri('bilbao',256));
%         %  ax.XTickLabel = [];
%         % ax.YTickLabel = [];
%         % Select Data

        % Produce plot and Color with final Lambda
% F2 = figure(5);
%         set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%         hold on
%         plot(t_dim,D(1:length(t)),"Color",cmap(index_color1,:),LineWidth=0.8,LineStyle="-")
%         ax2 = gca;
% 
%         %ax.XLim = [0,10];
%         %ax.YLim = [-1.5,1.5];
%         %ax2.XLabel.String = '$t [\mathrm{Myrs}]$';
%         %ax2.YLabel.String = '${H} [\mathrm{km}]$';
%         ax2.XLabel.Interpreter = 'latex';
%         ax2.YLabel.Interpreter = 'latex';
%         ax2.FontSize = font_axes;
%         ax2.TickLabelInterpreter = 'latex';
%         ax2.XMinorTick = 'on';
%         ax2.YMinorTick = 'on';
% 
%         ax2.Box  = 'on';
%         ax2.XGrid = 'on';
%         ax2.YGrid = 'on';
%         ax2.XColor = [0 0 0 ];
%         ax2.YColor = [0 0 0 ];
%         ax2.LineWidth = 1.2;
% 
%         cmap = colormap(crameri('bilbao',256));
%         %  ax.XTickLabel = [];
%         % ax.YTickLabel = [];
%         % Select Data

        % Produce plot and Color with final Lambda



        % Thickness
%         F3 = figure(6);
%         set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%         hold on
%         cmap = colormap(crameri('oslo',256));
% 
%         plot(t,D,"Color",cmap(index_color2,:),LineWidth=0.8,LineStyle="-")
%         ax2 = gca;
%         ax2.XLim = [0,10];
%         ax2.YLim = [0.1,1.1];
%         % ax2.XLabel.String = '$t^{\dagger}$';
%         %ax2.YLabel.String = '${D}^{\dagger, 2D}$';
%         ax2.XLabel.Interpreter = 'latex';
%         ax2.YLabel.Interpreter = 'latex';
%         ax2.FontSize = font_axes;
%         ax2.TickLabelInterpreter = 'latex';
%         ax2.Box  = 'on';
%         ax2.XGrid = 'on';
%         ax2.YGrid = 'on';
%         ax2.XColor = [0 0 0 ];
%         ax2.YColor = [0 0 0 ];
%         ax2.LineWidth = LineWidth;
%         %  ax2.XTickLabel = [];
%         % ax2.YTickLabel = [];
% 
%         F4 = figure(7);
%         set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%         hold on
%         cmap = colormap(crameri('oslo',256));
% 
%         plot(t,D_0D,"Color",cmap(index_color2,:),LineWidth=0.8,LineStyle="-")
%         ax2 = gca;
%         ax2.XLim = [0,10];
%         ax2.YLim = [0.1,1.1];
%         % ax2.XLabel.String = '$t^{\dagger}$';
%         % ax2.YLabel.String = '${D}^{\dagger, 0D}$';
%         ax2.XLabel.Interpreter = 'latex';
%         ax2.YLabel.Interpreter = 'latex';
%         ax2.FontSize = font_axes;
%         ax2.TickLabelInterpreter = 'latex';
%         ax2.Box  = 'on';
%         ax2.XGrid = 'on';
%         ax2.YGrid = 'on';
%         ax2.XColor = [0 0 0 ];
%         ax2.YColor = [0 0 0 ];
%         ax2.LineWidth = LineWidth;
%         %  ax2.XTickLabel = [];
%         % ax2.YTickLabel = [];
%         %TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
%         %color with residuum
% 
%         % Tau
%         F5 = figure(8);
%         set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%         hold on
%         cmap = colormap(crameri('oslo',256));
% 
%         plot(t,tau(1:length(t)),"Color",cmap(index_color2,:),LineWidth=0.8,LineStyle="-")
%         ax2 = gca;
%         ax2.XLim = [0,10];
%         %ax.YLim = [-1.5,1.5];
%         % ax2.XLabel.String = '$t^{\dagger}$';
%         %ax2.YLabel.String = '$\tau_{eff}^{\dagger, 2D}$';
%         ax2.XLabel.Interpreter = 'latex';
%         ax2.YLabel.Interpreter = 'latex';
%         ax2.FontSize = font_axes;
%         ax2.TickLabelInterpreter = 'latex';
%         ax2.Box  = 'on';
%         ax2.XGrid = 'on';
%         ax2.YGrid = 'on';
%         ax2.XColor = [0 0 0 ];
%         ax2.YColor = [0 0 0 ];
%         ax2.LineWidth = LineWidth;
%         %  ax2.XTickLabel = [];
%         % ax2.YTickLabel = [];
% 
%         F6 = figure(8);
%         set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%         hold on
%         cmap = colormap(crameri('oslo',256));
% 
%         plot(t0D,tau0D(3,:),"Color",cmap(index_color2,:),LineWidth=0.8,LineStyle="-")
%         ax2 = gca;
%         ax2.XLim = [0,10];
%         %ax.YLim = [-1.5,1.5];
%         %ax2.XLabel.String = '$t^{\dagger}$';
%         %ax2.YLabel.String = '$\tau_{eff}^{\dagger, 0D}$';
%         ax2.XLabel.Interpreter = 'latex';
%         ax2.YLabel.Interpreter = 'latex';
%         ax2.FontSize = font_axes;
%         ax2.TickLabelInterpreter = 'latex';
%         ax2.Box  = 'on';
%         ax2.XGrid = 'on';
%         ax2.YGrid = 'on';
%         ax2.XColor = [0 0 0 ];
%         ax2.YColor = [0 0 0 ];
%         ax2.LineWidth = LineWidth;
%         %  ax2.XTickLabel = [];
%         % ax2.YTickLabel = [];
%         % color with residuum

        % Topography
        F100 = figure(100);
        cmap = colormap(lipari);

        set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

        hold on
        plot(t_M,dU./(abs(min(dU))),"Color",cmap(index_color1,:),LineWidth=0.8,LineStyle="-")
        ax = gca;

        ax.XLim = [0,10];
        ax.YLim = [-0.1,0.1];
        %ax.XLabel.String = '$t^{\dagger}$';
        %ax.YLabel.String = '${H} [\mathrm{km}]$';
        ax.XLabel.Interpreter = 'latex';
        ax.YLabel.Interpreter = 'latex';
        ax.FontSize = font_axes;
        ax.TickLabelInterpreter = 'latex';
        ax.Box  = 'on';
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.XColor = [0 0 0 ];
        ax.YColor = [0 0 0 ];
        ax.LineWidth = 1.0;
        % ax.XTickLabel =[];
        %ax.YTickLabel = [];
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';



        F101 = figure(101);
        cmap = colormap(lipari);

        set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])

        hold on
        plot(T_Mean,dU,"Color",cmap(index_color1,:),LineWidth=0.8,LineStyle="-")
        ax = gca;

        % ax.XLim = [0,10];
        ax.YLim = [-0.1,0.2];
        %ax.XLabel.String = '$t^{\dagger}$';
        %ax.YLabel.String = '${H} [\mathrm{km}]$';
        ax.XLabel.Interpreter = 'latex';
        ax.YLabel.Interpreter = 'latex';
        ax.FontSize = font_axes;
        ax.TickLabelInterpreter = 'latex';
        ax.Box  = 'on';
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.XColor = [0 0 0];%[0 102 204 ]./255;
        ax.YColor = [0 0 0];%[0 102 204 ]./255;
        ax.LineWidth = 1.4;
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        %ax.XTickLabel =[];
        %ax.YTickLabel = [];

    end
end
F1;
figure_name ='Topo_ND';

pt=fullfile(pt_save,figure_name);
print(F1,pt,'-dpng')


% F2;
% figure_name ='Topo_D';
% 
% pt=fullfile(pt_save,figure_name);
% print(F2,pt,'-dpng')

% F3;
% figure_name ='D_2D';
% 
% pt=fullfile(pt_save,figure_name);
% print(F3,pt,'-dpng')
% F4;
% figure_name ='D_0D';

% pt=fullfile(pt_save,figure_name);
% print(F4,pt,'-dpng')
% F5;
% figure_name ='tau_2D';
% 
% pt=fullfile(pt_save,figure_name);
% print(F5,pt,'-dpng')
% F6;
% figure_name ='tau_0D';
% 
% pt=fullfile(pt_save,figure_name);
% print(F6,pt,'-dpng')
% 
% F7;
% figure_name ='tau_0L';
% pt=fullfile(pt_save,figure_name);
% print(F7,pt,'-dpng')
% 
% F100;
% 
% figure_name ='dU_1';
% pt=fullfile(pt_save,figure_name);
% print(F100,pt,'-dpng');

F101;

figure_name ='dU_2';
pt=fullfile(pt_save,figure_name);
print(F101,pt,'-dpng');

figure(100) %ColorBar fit
clf;

cmap = colormap(crameri('oslo',256));
ax = axes;
c = colorbar(ax);
ax.Visible = 'off';
Ticks=([0.0:5.0:c_l2(end)]);
TickLabels=arrayfun(@(x) sprintf('${%.2f}$',x),Ticks,'UniformOutput',false);

Ticks=([0.0:5.0:c_l2(end)]);
% Produce the tick
%TickLabels=arrayfun(@(x) sprintf('${%.2f}$',x),Ticks,'UniformOutput',false);

%TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
c.Ticks= Ticks;
c.TickLabelInterpreter = 'latex';
c.TickLabels=TickLabels;
%c2.Location = 'southoutside';
c.FontSize = 13;
%c.Label.String = 'F, [\%]';
c.Color    = [0,0,0];
c.Label.Interpreter = 'latex';

caxis(c_l2)
pt=fullfile(pt_save,'cbar_fit');
print(pt,'-dpng')

figure(100) %Colorbar Lambda
clf;
cmap = colormap(lipari);
ax = axes;
c = colorbar(ax);
ax.Visible = 'off';
Ticks=round([min(c_l1):max(c_l1)]);
% Produce the tick
%TickLabels=arrayfun(@(x) sprintf('$10^{%d}$',x),Ticks,'UniformOutput',false);
%TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
c.Ticks= Ticks;
%c.TickLabelInterpreter = 'latex';
%c.TickLabels=TickLabels;
%c.Location = 'southoutside';
c.LineWidth = 1.2;
c.FontSize = 13;
%c.Label.String = '$\frac{\Lambda}{1+\xi^{UM}}$';
c.Label.Interpreter = 'latex';
caxis(c_l1);
pt=fullfile(pt_save,'cbar_topo');
print(pt,'-dpng')
close all;
end

function [index_cmap] = find_color(V,clim,islog)
z_min = min(clim);
z_max = max(clim);
if islog == 1
    V = log10(V);
end
V = (V-z_min)/(z_max-z_min);
if V<0
    V=0;
elseif V>1
    V=1;
end
V=round(1+V*(256-1));%round to nearest index
index_cmap=V;
end


function [II]= plot_2Dmaps_profiles(S,pt_save,nlm)
close all;
clf;
FIT = S.FIT;
TB  = S.TB;
font_axes = 16;
font_legend = 14;
font_text   = 5;
size_shit = [12,13.5];
LineWidth = 1.2;
marker_size = 10;
fnames = fieldnames(TB);
z = -1000:1050/512:50;
z = z(z<20);
ind_z = find(z>-200);
flength = length(fnames);

for ktest = 1:flength
    if TB.(fnames{ktest}).P_Var.failed == 0
        t = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).ID.n;
        D_m = TB.(fnames{ktest}).D_matrix;
        D_m(D_m==-Inf)=nan;
        for i=1:length(t)
            a = D_m(i,:);
            D_mv(i,:) = movmean(a,10);
            b = D_mv(i,:);
            ind=find(b(ind_z(1):ind_z(end))==nanmin(b(ind_z(1):ind_z(end))),1);
            x(i) = z(ind_z(ind(1)));
            DD(i)=b(ind_z(ind(1)));
        end
        if strcmp(fnames{ktest},'T22')
        bla=0; 
        end
        figure(1)
        clf;
        ax = gca;
        pcolor(t,z,D_m(1:length(t),:)'./80);shading flat; colorbar;
        colorbar;
        caxis([0.1,1.0])
        colormap(crameri('nuuk',27))
        hold on
        ylim([-300,-100])
        ax.Box = 'on';
        ax.LineWidth = 1.2;
        ax.XColor = [0 0 0 ];
        ax.YColor = [0 0 0 ]; 
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.Layer = 'top';
        figure_name = fnames{ktest};
        pt_fig = fullfile(pt_save,figure_name);
        print(pt_fig,'-dpng','-r600');
        bla =0;

        figure(2)
        clf;
        ax = gca; 
        hold on

        plot(D_mv./80,z,'k')
        plot(DD./80,x,'r','LineWidth',1.2)
        yline(TB.(fnames{ktest}).Dp2-100,Color='b',LineStyle=':',LineWidth=1.2)
        yline(-100-TB.(fnames{ktest}).f(1).*TB.(fnames{ktest}).P_Var.L0./1e3,Color='r',LineStyle='--',LineWidth=1.2);
        hold on
        yline(TB.(fnames{ktest}).Dp2-100,Color='b',LineStyle=':',LineWidth=1.2)
        yline(-100-TB.(fnames{ktest}).f(1).*TB.(fnames{ktest}).P_Var.L0./1e3,Color='r',LineStyle='--',LineWidth=1.2)
        xlim([0.01,0.9])
        ylim([-300,-100])
        ax.Box = 'on';
        ax.LineWidth = 1.2;
        ax.XColor = [0 0 0 ];
        ax.YColor = [0 0 0 ]; 
        ax.XGrid = 'on';
        ax.YGrid = 'on'; 
        figure_name = strcat(figure_name,'profile');
        pt_fig = fullfile(pt_save,figure_name);
        print(pt_fig,'-dpng','-r600');

        figure(3)
        hold on 
        buf1 = -100-TB.(fnames{ktest}).f(1).*TB.(fnames{ktest}).P_Var.L0./1e3;
        buf2 = TB.(fnames{ktest}).Dp2-100; 
        I = (x-buf1)./(mean(x(end-1))-buf1); 
        I = movmean(I,4);
        plot(t,log10(abs(I)))
        II(ktest)=(buf1-mean(x(end-1)))./(mean(x(end-1)));

        x = [];
        DD= [];
        D_mv = [];
        I = [];

    end
end


end
