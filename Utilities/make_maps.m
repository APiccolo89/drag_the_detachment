function make_maps(x,y,z,cmap,logx,logy,label1,label2,label3,filename,ptsave,nfigure)

figure(nfigure)
ax = gca; 
p1 = contourf(x,y,z,10);
colormap(cmap)
shading interp; 
ax.XScale = logx;
ax.YScale = logy;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
ax.XLabel.String = label1;
ax.YLabel.String = label2;
ax.LineWidth = 1.2; 
ax.Box = 'on';
ax.Layer = 'top';
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.TickDir = 'both';
cb = colorbar(ax);
cb.Label.Interpreter = 'latex';
cb.Label.String     = label3; 
cb.Box = 'on';
cb.TickLabelInterpreter = 'latex';

end

