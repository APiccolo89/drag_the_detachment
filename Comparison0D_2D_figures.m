clear all;
close all;
addpath('Utilities/')
addpath('Utilities/Plot_Class/')
addpath('Utilities/ScientificColourMaps8/')
% Comparison 0D - 2D
load('Data_Base_Fit.mat')

path_save = '../../Manuscript_Revision';


%plot_comparison(NLN.TB.T1,NLN2.TB.T1,path_save,'Comparison_T1',[0,3.0])
%plot_comparison(NLN.TB.T15,NLN2.TB.T15,path_save,'Comparison_T15',[0,9.0])
% Plot data Linear Experiment database
% create path where to save if does not exit
nlm = Problem_type;
nlm.Linear = 1;

 pt = fullfile(path_save,'Maps_Linear');
 if ~isdir(pt)
     mkdir(pt)
 end
% %%
% % Create figure
plot_the_figures(LN_fetches1,pt,nlm)
%[linear]=plot_2Dmaps_profiles(LN_fetches1,pt,nlm);
% %nlm.Linear = 0;
% %[Nlinear]=plot_2Dmaps_profiles(NLN_fetches2,pt,nlm);
% pt = fullfile(path_save,'Plots_2Fetches');
 %plot_the_figures(NLN_fetches2,pt,nlm)
% disp('===========================================================')
% disp('Non Linear 2 fetches')
% disp('===========================================================')
% 
[Nlinear]=plot_2Dmaps_profiles(NLN_fetches2,pt,nlm);
% pt = fullfile(path_save,'Plots_1Fetches');
% plot_the_figures(NLN_fetches1,pt,nlm)
% %disp('===========================================================')
% %disp('Non Linear 1 fetches')
% %disp('===========================================================')
% 
% 
% % % Plot data Non Linear Experiment with one fetches
% pt = fullfile(path_save,'Plots_linear');
% nlm.Linear=1;
% plot_the_figures(LN_fetches1,pt,nlm)
% 
% disp('===========================================================')
% disp(' Linear 1 fetches')
% disp('===========================================================')
pt1 = fullfile(path_save,'Maps_Linear');
pt2 = fullfile(path_save,'Maps_NLinear');
pt3 = fullfile(path_save,'Maps_NLinear2');


[linear]=plot_2Dmaps_profiles(LN_fetches1,pt1,nlm);
[Nlinear2]=plot_2Dmaps_profiles(NLN_fetches2,pt2,nlm);
[Nlinear1]=plot_2Dmaps_profiles(NLN_fetches1,pt3,nlm);

%%
% box figure
pt = path_save;
figure(2)

clf;
ax = gca; 
Data = [ones(length(linear),1)'];
Data2 = [ones(length(Nlinear1),1)'.*2];
Data3 = [ones(length(Nlinear2),1)'.*3];
hold on; 
boxchart(Data,linear,"BoxFaceColor",[ 230, 126, 34]./255,'Notch','on')
boxchart(Data2,Nlinear1,"BoxFaceColor",[ 33, 97, 140]./255,'Notch','on')
boxchart(Data3,Nlinear2,"BoxFaceColor",[ 89, 218, 120]./255,'Notch','on')

ax.Box = 'on';
ax.LineWidth = 1.2; 
%ax.XGrid = 'on';
%ax.YGrid = 'on';
ax.XColor = [0 0 0 ];
ax.YColor = [0 0 0 ]; 
figure_name = 'box_plot.png';
pt = fullfile(path_save,figure_name);
print(pt,'-dpng','-r600')

%%
plot_comparison(NLN_fetches1.TB.T3,NLN_fetches2.TB.T1,path_save,'Comparison_T3',[0,3.0])
plot_comparison(NLN_fetches1.TB.T15,NLN_fetches2.TB.T15,path_save,'Comparison_T15',[0,9.0])
%%
function plot_comparison(Tests,Test2fetches,pt_save,name_figure,tlim)


clf;
close all;
Lambda_0 = Tests.ID.Lambda./(1+Tests.ID.Df_UM.*(Tests.ID.tau_mc).^(3.5-1));

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
%ax2.XGrid = 'on';
%ax2.YGrid = 'on';
ax2.XColor = [0 0 0 ];
ax2.YColor = [0 0 0 ];
ax2.LineWidth = 1.2;
l=legend('$2D$','$0D(\omega_{depth})$','$0D(\omega_S,\omega_{depth})$','Analytical solution');
l.Interpreter = 'latex';
l.FontSize = font_legend;
l.Location = 'northeast';
str = ['$\mathbf{log_{10}\left(\Lambda_0\right)} = ',num2str(log10(Lambda_0),3) ,'$'];
text = annotation('textbox',[0.125,0.85,0.2,0.05],'string',str);
text.FontSize = 10;
text.Interpreter = 'latex';
text.LineStyle = "none";


pt=fullfile(pt_save,name_figure);
print(pt,'-dpng')

end

