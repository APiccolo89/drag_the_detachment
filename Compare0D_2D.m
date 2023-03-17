%=========================================================================%

clear all;
close all;
addpath Adimensional\
addpath Dimensional\
addpath \Users\'Andrea Piccolo'\Dropbox\freezeColors-master\
addpath Utilities\
addpath D0_D2_comparison_function\
clear all
% Common Data
D0=80e3;
n=3.5;
Df_S=10;
nlm = Problem_type.NonLinear;
Df_UM = nan;
DB_path = '../Test_DB2.hdf5';
pt_save = '../../Results/PRESENTATION_Comparison2D_NL';
if not(isdir(pt_save))
     mkdir(pt_save);
end
l = h5info(DB_path,'/Viscous/HR_DIS');
%l2 =h5info(DB_path,'/Viscous/LR_NLM');
TestsA  = {l.Groups.Name};
%TestsB  = {l2.Groups.Name};
HR = 1:length(TestsA);
%LR = -(1:length(TestsB));
Res = [HR]; 
Tests   = [TestsA]; 
%function Optimize Data Base 

[TB,FIT] = perform_optimization_DataBase(Tests,n,D0,Df_S,nlm,Df_UM,DB_path,pt_save,Res);
%%
%clf; 
%close all; 

figure(5)
ax = gca; 
scatter(FIT.Lambda,-FIT.Dp./2,20,[204 204 255]./255,'filled','o','MarkerEdgeColor','k');
hold on
scatter(FIT.Lambda,FIT.fetch(1,:),20,[204 255 204]./255,'filled','d','MarkerEdgeColor','k')
ax.XScale = 'log';
ax.XColor = [0,0,0];
ax.XLabel.String = "$log\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$"; 
ax.XLabel.Interpreter = 'latex';
ax.YColor = [0,0,0];
ax.YLabel.String = "$\omega_{\mathrm{depth}}$"; 
ax.YLabel.Interpreter = 'latex';
ax.XGrid = 'on';
ax.XMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YGrid = 'on';
ax.YMinorTick = 'on';
ax.YMinorGrid = 'on';
ax.Box = 'on';
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.0; 
ax.FontSize = 16; 
l = legend('2D simulations','0D fitting');
l.Interpreter = 'latex'; 
figure_name ='Depth_Coefficient'; 

pt=fullfile(pt_save,figure_name);
print(pt,'-dpng')
figure(7)
ax = gca; 
scatter(FIT.Lambda,FIT.fetch(2,:),20,[204 204 255]./255,'filled','o','MarkerEdgeColor','k');
hold on
ax.XScale = 'log';
ax.YScale = 'log';

ax.XColor = [0,0,0];
ax.YLim = [0,100];
ax.XLabel.String = "$log\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$"; 
ax.XLabel.Interpreter = 'latex';
ax.YColor = [0,0,0];
ax.YLabel.String = "$\omega_{\mathrm{S}}$"; 
ax.YLabel.Interpreter = 'latex';
ax.XGrid = 'on';
ax.XMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YGrid = 'on';
ax.YMinorTick = 'on';
ax.YMinorGrid = 'on';
ax.Box = 'on';
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.0; 
ax.FontSize = 16; 
figure_name ='S_Coefficient'; 

pt=fullfile(pt_save,figure_name);
print(pt,'-dpng')






figure(8)
ax = gca; 
scatter(FIT.Lambda,FIT.fitting_p,20,[204 204 255]./255,'filled','o','MarkerEdgeColor','k');
ax.XScale = 'log';
ax.XColor = [0,0,0];
ax.XLabel.String = "$log\left(\frac{\Lambda}{1+\xi^{\mathrm{UM}}}\right)$"; 
ax.XLabel.Interpreter = 'latex';
ax.YColor = [0,0,0];
ax.YLabel.String = "$F$"; 
ax.YLabel.Interpreter = 'latex';
ax.XGrid = 'on';
ax.XMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YGrid = 'on';
ax.YMinorTick = 'on';
ax.YMinorGrid = 'on';
ax.Box = 'on';
ax.TickLabelInterpreter = 'latex';
ax.LineWidth = 1.0; 
ax.FontSize = 16; 
%l = legend('2D simulations','0D fitting');
%l.Interpreter = 'latex'; 
figure_name ='Fitting'; 

pt=fullfile(pt_save,figure_name);
print(pt,'-dpng')


%% 1D line plot 
fnames = fieldnames(TB);
flength = length(fnames); 

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

    if ~isempty(D_0D)
    c_l1 = [-8,0];% lim color Lambda
    c_l2 = [0,2];% lim color Res
    R    = TB.(fnames{ktest}).res.*100;
    L    = TB.(fnames{ktest}).ID.Lambda; 
    if nlm.islinear == 0
        L = L./(TB.(fnames{ktest}).ID.Df_UM);
    end
    index_color1 = find_color(L,c_l1,1);
    index_color2 = find_color(R,c_l2,1);

    cmap = colormap(crameri('bilbao',256));
    % Topography
    F1 = figure(9);
    hold on
    plot(t,topo(1:length(t)),"Color",cmap(index_color1,:),LineWidth=0.8,LineStyle="-")
    ax = gca;
    
    ax.XLim = [0,10];
    %ax.YLim = [-1.5,1.5];
    ax.XLabel.String = '$t^{\dagger}$';
    ax.YLabel.String = '${H} [\mathrm{km}]$';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.FontSize = 15;
    ax.TickLabelInterpreter = 'latex';
    ax.Box  = 'on';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.XColor = [0 0 0 ];
    ax.YColor = [0 0 0 ];
    ax.LineWidth = 1.0;
    cmap = colormap(crameri('bilbao',256));
    c = colorbar;
    Ticks=round([min(c_l1):max(c_l1)]);
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('$10^{%d}$',x),Ticks,'UniformOutput',false);
    %TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
    c.Ticks= Ticks;
    c.TickLabelInterpreter = 'latex';
    c.TickLabels=TickLabels;
    %c.Location = 'eastoutside';
    c.FontSize = 12;
    c.Label.String = '$\frac{\Lambda}{1+\xi^{\mathrm{UM}}}$';
    c.Label.Interpreter = 'latex';
    caxis(c_l1);
    

    F2 = figure(10);
    hold on
    plot(t_dim,topo(1:length(t)),"Color",cmap(index_color1,:),LineWidth=0.8,LineStyle="-")
    ax2 = gca;
    %ax.XLim = [0,10];
    %ax.YLim = [-1.5,1.5];
    ax2.XLabel.String = '$t [\mathrm{Myrs}]$';
    ax2.YLabel.String = '${H} [\mathrm{km}]$';
    ax2.XLabel.Interpreter = 'latex';
    ax2.YLabel.Interpreter = 'latex';
    ax2.FontSize = 15;
    ax2.TickLabelInterpreter = 'latex';
    ax2.Box  = 'on';
    ax2.XGrid = 'on';
    ax2.YGrid = 'on';
    ax2.XColor = [0 0 0 ];
    ax2.YColor = [0 0 0 ];
    ax2.LineWidth = 1.0;
    cmap = colormap(crameri('bilbao',256));
    c = colorbar;freezeColors;
    Ticks=round([min(c_l1):max(c_l1)]);
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('$10^{%d}$',x),Ticks,'UniformOutput',false);
    %TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
    c.Ticks= Ticks;
    c.TickLabelInterpreter = 'latex';
    c.TickLabels=TickLabels;
    %c.Location = 'eastoutside';
    c.FontSize = 12;
    c.Label.String = '$\frac{\Lambda}{1+\xi^{\mathrm{UM}}}$';
    c.Label.Interpreter = 'latex';
    caxis(c_l1);
    % Select Data 

    % Produce plot and Color with final Lambda



    % Thickness
    F3 = figure(11);
      hold on
    cmap = colormap(crameri('oslo',256));

    plot(t,D,"Color",cmap(index_color2,:),LineWidth=0.8,LineStyle="-")
    ax2 = gca;
    ax2.XLim = [0,10];
    ax2.YLim = [0.1,1.1];
    ax2.XLabel.String = '$t^{\dagger}$';
    ax2.YLabel.String = '${D}^{\dagger, 2D}$';
    ax2.XLabel.Interpreter = 'latex';
    ax2.YLabel.Interpreter = 'latex';
    ax2.FontSize = 15;
    ax2.TickLabelInterpreter = 'latex';
    ax2.Box  = 'on';
    ax2.XGrid = 'on';
    ax2.YGrid = 'on';
    ax2.XColor = [0 0 0 ];
    ax2.YColor = [0 0 0 ];
    ax2.LineWidth = 1.0;
    c2 = colorbar;
    Ticks=round([min(c_l2):1:max(c_l2)]);
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('${%d}$',x),Ticks,'UniformOutput',false);
    %TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
    c2.Ticks= Ticks;
    c2.TickLabelInterpreter = 'latex';
    c2.TickLabels=TickLabels;
    %c.Location = 'eastoutside';
    c2.FontSize = 12;
    ax2.Color=[.90,.90,.90];
    c2.Label.String = 'F, [%]';
   
    c2.Color    = [0,0,0];
    caxis(c_l2)


    F4 = figure(12);
      hold on
    cmap = colormap(crameri('oslo',256));

    plot(t,D_0D,"Color",cmap(index_color2,:),LineWidth=0.8,LineStyle="-")
    ax2 = gca;
    ax2.XLim = [0,10];
    ax2.YLim = [0.1,1.1];
    ax2.XLabel.String = '$t^{\dagger}$';
    ax2.YLabel.String = '${D}^{\dagger, 0D}$';
    ax2.XLabel.Interpreter = 'latex';
    ax2.YLabel.Interpreter = 'latex';
    ax2.FontSize = 15;
    ax2.TickLabelInterpreter = 'latex';
    ax2.Box  = 'on';
    ax2.XGrid = 'on';
    ax2.YGrid = 'on';
    ax2.XColor = [0 0 0 ];
    ax2.YColor = [0 0 0 ];
    ax2.LineWidth = 1.0;
    c2 = colorbar;
    Ticks=round([min(c_l2):1:max(c_l2)]);
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('${%d}$',x),Ticks,'UniformOutput',false);
    %TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
    c2.Ticks= Ticks;
    c2.TickLabelInterpreter = 'latex';
    c2.TickLabels=TickLabels;
    %c.Location = 'eastoutside';
    c2.FontSize = 12;
    ax2.Color=[.90,.90,.90];
   c2.Label.String = 'F, [%]';
  %  c2.Label.Interpreter = 'latex';
    c2.Color    = [0,0,0];
    caxis(c_l2)


    %color with residuum 

    % Tau 
    F5 = figure(13);
          hold on
    cmap = colormap(crameri('oslo',256));

    plot(t,tau(1:length(t)),"Color",cmap(index_color2,:),LineWidth=0.8,LineStyle="-")
    ax2 = gca;
    ax2.XLim = [0,10];
    %ax.YLim = [-1.5,1.5];
    ax2.XLabel.String = '$t^{\dagger}$';
    ax2.YLabel.String = '$\tau_{eff}^{\dagger, 2D}$';
    ax2.XLabel.Interpreter = 'latex';
    ax2.YLabel.Interpreter = 'latex';
    ax2.FontSize = 15;
    ax2.TickLabelInterpreter = 'latex';
    ax2.Box  = 'on';
    ax2.XGrid = 'on';
    ax2.YGrid = 'on';
    ax2.XColor = [0 0 0 ];
    ax2.YColor = [0 0 0 ];
    ax2.LineWidth = 1.0;
    c2 = colorbar;
    Ticks=round([min(c_l2):1:max(c_l2)]);
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('${%d}$',x),Ticks,'UniformOutput',false);
    %TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
    c2.Ticks= Ticks;
    c2.TickLabelInterpreter = 'latex';
    c2.TickLabels=TickLabels;
    %c.Location = 'eastoutside';
    c2.FontSize = 12;
    ax2.Color=[.90,.90,.90];
   c2.Label.String = 'F, [%]';
  %  c2.Label.Interpreter = 'latex';
    c2.Color    = [0,0,0];
    caxis(c_l2)

    F6 = figure(14);
              hold on
    cmap = colormap(crameri('oslo',256));

    plot(t0D,tau0D(3,:),"Color",cmap(index_color2,:),LineWidth=0.8,LineStyle="-")
    ax2 = gca;
    ax2.XLim = [0,10];
    %ax.YLim = [-1.5,1.5];
    ax2.XLabel.String = '$t^{\dagger}$';
    ax2.YLabel.String = '$\tau_{eff}^{\dagger, 0D}$';
    ax2.XLabel.Interpreter = 'latex';
    ax2.YLabel.Interpreter = 'latex';
    ax2.FontSize = 15;
    ax2.TickLabelInterpreter = 'latex';
    ax2.Box  = 'on';
    ax2.XGrid = 'on';
    ax2.YGrid = 'on';
    ax2.XColor = [0 0 0 ];
    ax2.YColor = [0 0 0 ];
    ax2.LineWidth = 1.0;
    c2 = colorbar;
    Ticks=round([min(c_l2):1:max(c_l2)]);
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('${%d}$',x),Ticks,'UniformOutput',false);
    %TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
    c2.Ticks= Ticks;
    c2.TickLabelInterpreter = 'latex';
    c2.TickLabels=TickLabels;
    %c.Location = 'eastoutside';
    c2.FontSize = 12;
    ax2.Color=[.90,.90,.90];
   c2.Label.String = 'F, [%]';
    c2.Color    = [0,0,0];
    caxis(c_l2)
    % color with residuum 
    end
end
F1;
figure_name ='Topo_ND'; 

pt=fullfile(pt_save,figure_name);
print(F1,pt,'-dpng')


F2;
figure_name ='Topo_D'; 

pt=fullfile(pt_save,figure_name);
print(F2,pt,'-dpng')

F3;
figure_name ='D_2D'; 

pt=fullfile(pt_save,figure_name);
print(F3,pt,'-dpng')
F4;
figure_name ='D_0D'; 

pt=fullfile(pt_save,figure_name);
print(F4,pt,'-dpng')
F5;
figure_name ='tau_2D'; 

pt=fullfile(pt_save,figure_name);
print(F5,pt,'-dpng')
F6;
figure_name ='tau_0D'; 

pt=fullfile(pt_save,figure_name);
print(F6,pt,'-dpng')


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

