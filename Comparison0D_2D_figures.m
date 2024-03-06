clear all;
close all;
addpath('Utilities/')
addpath('Utilities/Plot_Class/')
addpath('Utilities/ScientificColourMaps8/')
% Comparison 0D - 2D
load('Data_Base_Fit_DEF.mat')

path_save = 'New_Manuscript_Figure';


plot_comparison(NLN_fetches2.TB.T_0_1,NLN_fetches2.TB.T_0_15,path_save,'Comparison_T1.png',[0,3.0])
%plot_comparison(NLN_fetches2.TB.T_0_15,path_save,'Comparison_T15',[0,9.0])
% Plot data Linear Experiment database
% create path where to save if does not exit
nlm = Problem_type;
nlm.Linear = 1;

 pt = fullfile(path_save,'Maps_Linear');
 if ~isdir(pt)
     mkdir(pt)
 end

%figure_5(NLN_fetches2,path_save)

%figure_6(NLN_fetches2,path_save)























%%
function plot_comparison(Tests1,Tests2,pt_save,name_figure,tlim)


clf;
close all;
Lambda_01 = Tests1.ID.Lambda./(1+Tests1.ID.Df_UM.*(Tests1.ID.tau_mc).^(3.5-1));
Lambda_02 = Tests2.ID.Lambda./(1+Tests2.ID.Df_UM.*(Tests2.ID.tau_mc).^(3.5-1));



F1 = figure('DefaultTextFontName', 'Comic Sans MS', 'DefaultAxesFontName', 'Comic Sans MS');

t_d = 0:0.001:1;
D   = (1-t_d).^(1/3.5);

set(gcf, 'Units','centimeters', 'Position', [0, 0, 7.5,13], 'PaperUnits', 'centimeters', 'PaperSize', [8, 12])
subplot(2,1,1)

plot(Tests1.t.*3.5,Tests1.D,"Color",'#9a2a3c',LineWidth=1.2,LineStyle="-")
hold on
%plot(Tests.t.*3.5,Tests.D0D2,'Color','#b2ec73',LineWidth=1.0,LineStyle='--')
%plot(Test2fetches.t.*3.5,Test2fetches.D0D2,'Color','#6e4e8e',LineWidth=1.0,LineStyle='--')
plot(t_d,D,'Color','k','LineStyle','-.','LineWidth',1.0)

ax2 = gca;
ax2.Position = [0.1875    0.5500    0.7750    0.4000];
ax2.XLim = [0,10];
ax2.YLim = [0.1,1.1];
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.Interpreter = 'latex';
ax2.FontUnits = 'centimeters';
ax2.XTick      = [0,2,4,6,8,10];
ax2.XTickLabel = {'0',[],'4',[],'8',[]};
ax2.YTick = [0.1,0.5,1.0];
ax2.YTickLabel =  {'0.1',[],'1.0'};
ax2.FontSize = 0.5;
%ax2.FontSize = 0.5;
ax2.TickLabelInterpreter = 'latex';
ax2.Box  = 'on';
%ax2.XGrid = 'on';
%ax2.YGrid = 'on';
ax2.XColor = [0 0 0 ];
ax2.YColor = [0 0 0 ];
ax2.LineWidth = 1.2;
ax2.XLabel.String = [];
ax2.XTickLabel = [];
ax2.YLabel.String = '$D^{\dagger}$';
ax2.YLabel.Interpreter = 'latex'; 
ax2.YLabel.Units = 'normalized';
ax2.YLabel.Position(2) = 0.5;
ax2.YLabel.Position(1) = -0.01;


l=legend('$2D$','Analytical solution');
l.Interpreter = 'latex';
l.FontSize = 8;
l.Location = 'northeast';
str = ['$\mathbf{log_{10}\left(\Lambda_0\right)} = ',num2str(log10(Lambda_01),3) ,'$'];
text = annotation('textbox',[0.125,0.85,0.2,0.05],'string',str);
text.FontSize = 8; 
text.Interpreter = 'latex';
text.LineStyle = "none";
text.Position(1) = 0.5;
text.Position(2) = 0.8;

subplot(2,1,2)
plot(Tests2.t.*3.5,Tests2.D,"Color",'#9a2a3c',LineWidth=1.2,LineStyle="-")
hold on
%plot(Tests.t.*3.5,Tests.D0D2,'Color','#b2ec73',LineWidth=1.0,LineStyle='--')
%plot(Test2fetches.t.*3.5,Test2fetches.D0D2,'Color','#6e4e8e',LineWidth=1.0,LineStyle='--')
plot(t_d,D,'Color','k','LineStyle','-.','LineWidth',1.0)

ax3 = gca;
ax3.Position = [    0.1875    0.1000    0.7750    0.4000];
ax3.XLim = [0,10.0];
ax3.YLim = [0.1,1.1];
ax3.XLabel.Interpreter = 'latex';
ax3.YLabel.Interpreter = 'latex';

ax3.TickLabelInterpreter = 'latex';
ax3.Box  = 'on';
%a32.XGrid = 'on';
%a32.YGrid = 'on';
ax3.XColor = [0 0 0 ];
ax3.YColor = [0 0 0 ];
ax3.LineWidth = 1.2;
ax3.XLabel.String = '$t^{\dagger}$';
ax3.XLabel.Interpreter = 'latex';
ax3.YLabel.String = '$D^{\dagger}$';
ax3.YLabel.Interpreter = 'latex'; 
l2=legend('$2D$','Analytical solution');
l2.Interpreter = 'latex';
l2.FontSize = 8;
l2.Location = 'northeast';
str = ['$\mathbf{log_{10}\left(\Lambda_0\right)} = ',num2str(log10(Lambda_02),3) ,'$'];
text2 = annotation('textbox',[0.125,0.85,0.2,0.05],'string',str);
text2.FontSize = 8;
text2.Interpreter = 'latex';
text2.LineStyle = "none";
text2.Position(1) = 0.5;
text2.Position(2) = 0.35;
ax3.FontUnits = 'centimeters';
ax3.XTick      = [0,2,4,6,8,10];
ax3.XTickLabel = {'0',[],'4',[],'8',[]};
ax3.YTick = [0.1,0.5,1.0];
ax3.YTickLabel =  {'0.1',[],'1.0'};
ax3.FontSize = 0.5;
ax3.YLabel.Units = 'normalized';
ax3.YLabel.Position(2) = 0.5;
ax3.YLabel.Position(1) = -0.01;
pt=fullfile(pt_save,name_figure);
exportgraphics(F1,pt,'Resolution',1000,'BackgroundColor','white')         
bla; 


end

