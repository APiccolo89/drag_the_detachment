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

