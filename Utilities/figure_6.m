function figure_6(Data_,ptsave)
% Create Data Set
%=========================================================================%
% Test_2D -> D,tau
% Test fitting -> D,tau
%=========================================================================%
FIT = Data_.FIT;

f=figure(1);
clf;

set(gcf, 'Units','centimeters', 'Position', [0, 0, 18,7], 'PaperUnits', 'centimeters', 'PaperSize', [15, 15])
subplot(1,3,1)
% Create the axis
ax0 = gca();
p0x  = 0.075;
p0y  = 0.15; 
sz_axh = 0.25;
sz_axv = 0.6;
dh     = sz_axh+0.06;
dv     = 0.00; 

ax0.Position = [p0x,p0y,sz_axh,sz_axv];
s1=scatter(ax0,FIT.Lambda,FIT.fetch(1,:),20,log10(FIT.xium),'filled','MarkerEdgeColor','k');
colormap(ax0,crameri('oslo'))
ax0.XScale = 'log';

ax1.XScale = 'log';
colormap(ax0,crameri('oslo'))
ax0.Box = 'on';
ax0.LineWidth = 1.2; 
ax0.FontUnits = 'centimeters';
ax0.FontSize = 0.5;
ax0.TickLabelInterpreter = 'latex';
ax0.XLabel.String = '$\Lambda_0$';
ax0.XLabel.Interpreter = 'latex';
ax0.XLim = [10^(-5),1];
ax0.XLabel.Interpreter = 'latex'; 
ax0.XTick = 10.^(-5:1:0);
ax0.XTickLabelRotation = 0;
ax0.XTickLabel = {'$10^{-5}$',[] ,[] ,[],'$10^{-1}$',[]};
ax0.YLabel.String = '$\omega_d$';
ax0.YLabel.Interpreter = 'latex';
ax0.YColor = 'k';
ax0.TickLabelInterpreter = 'latex';
ax0.YTick = [0,0.1,0.2,0.3,0.4];
ax0.YTickLabel = {'0.0',[],'0.2',[],'0.4'};
ax0.XLabel.Units = 'normalized';
ax0.YLabel.Units = 'normalized';
ax0.YLabel.Position(2) = 0.25;
ax0.YLabel.Position(1) = -0.01;
ax0.XLabel.Position(2) = -0.05;

subplot(1,3,2)
% Create the axis
ax1 = gca();
ax1.Position = [p0x+dh,p0y,sz_axh,sz_axv];
s2=scatter(ax1,FIT.Lambda,FIT.fetch(2,:),20,log10(FIT.xium),'filled','MarkerEdgeColor','k');
ax1.XScale = 'log';
colormap(ax1,crameri('oslo'))
ax1.Box = 'on';
ax1.LineWidth = 1.2; 
ax1.FontUnits = 'centimeters';
ax1.FontSize = 0.5;
ax1.XLabel.String = '$\Lambda_0$';
ax1.XLabel.Interpreter = 'latex';
ax1.XLim = [10^(-5),1];
ax1.XTick = 10.^(-5:1:0);
ax1.XTickLabelRotation = 0;
ax1.XTickLabel = {'$10^{-5}$',[] ,[] ,[],'$10^{-1}$',[]};
ax1.YLabel.String = '$\omega_s$';
ax1.YLabel.Interpreter = 'latex';
ax1.YColor = 'k';
ax1.TickLabelInterpreter = 'latex';
ax1.YTick = [0,1,2,3,4];
ax1.YTickLabel = {'0',[],'2',[],'4'};
ax1.XLabel.Units = 'normalized';
ax1.XLabel.Position(2) = -0.05;
ax1.YLabel.Units = 'normalized';
ax1.YLabel.Position(2) = 0.25;
ax1.YLabel.Position(1) = -0.01;


subplot(1,3,3)
% Create the axis
ax3 = gca();
ax3.Position = [p0x+2.*dh,p0y,sz_axh,sz_axv];
s3=scatter(ax3,FIT.Lambda,FIT.res_gf,20,log10(FIT.xium),'filled','MarkerEdgeColor','k');
ax3.XScale = 'log';
colormap(ax3,crameri('oslo'))
ax3.Box = 'on';
ax3.LineWidth = 1.2; 
ax3.FontUnits = 'centimeters';
ax3.FontSize = 0.5;
ax3.TickLabelInterpreter = 'latex';
ax3.XLabel.String = '$\Lambda_0$';
ax3.XLabel.Interpreter = 'latex';
ax3.XLim = [10^(-5),1];
ax3.XTick = 10.^(-5:1:0);
ax3.XTickLabelRotation = 0;
ax3.XTickLabel = {'$10^{-5}$',[] ,[] ,[],'$10^{-1}$',[]};
ax3.YLabel.String = 'Fit, $\%$';
ax3.YLabel.Interpreter = 'latex';
ax3.YColor = 'k';
ax3.TickLabelInterpreter = 'latex';
ax3.YTick = [0.0,0.05,0.1,0.15];
ax3.XLabel.Units = 'normalized';
ax3.XLabel.Position(2) = -0.05;
ax3.YTickLabel = {'0',[],'10',[]};
ax3.YLabel.Units = 'normalized';
ax3.YLabel.Position(2) = 0.3;
ax3.YLabel.Position(1) = -0.01;

cbar = colorbar(ax1,'northoutside');
cbar.Box = 'on';
cbar.TickLabelInterpreter = 'latex';
cbar.Label.String = '$log_{10}\left(\xi_M\right)$';
cbar.Label.Interpreter = 'latex';
cbar.Position = [p0x+sz_axh*0.5,0.85,2*dh,0.05];
cbar.TickLabelInterpreter = 'latex';
cbar.Label.Position(2)=-1.7;

t0 = annotation("textbox");
t0.Position = [p0x-0.01,p0y+dv+sz_axv-0.05,0.05,0.05];
t0.Interpreter = 'latex';
t0.String = '$[a]$';
t0.FontSize = 12;
t0.FontWeight = 'bold';

t0.LineStyle = 'none';

t1 = annotation("textbox");
t1.Position = [p0x+dh-0.01,p0y+dv+sz_axv-0.05,0.05,0.05];
t1.Interpreter = 'latex';
t1.String = '$[b]$';
t1.FontSize = 12;
t1.FontWeight = 'bold';
t1.LineStyle = 'none';


t2= annotation("textbox");
t2.Position = [p0x+2*dh-0.01,p0y+dv+sz_axv-0.05,0.05,0.05];
t2.Interpreter = 'latex';
t2.String = '$[c]$';
t2.FontSize = 12;
t2.FontWeight = 'bold';
t2.LineStyle = 'none';


figure_name = 'figure_6.png';
filename = fullfile(ptsave,figure_name);



exportgraphics(f,filename,'Resolution',1000,'BackgroundColor','white')         


end
