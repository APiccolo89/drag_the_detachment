function figure_5(Data_,ptsave)
% Create Data Set
%=========================================================================%
% Test_2D -> D,tau
% Test fitting -> D,tau
%=========================================================================%
TB = Data_.TB;
FIT = Data_.FIT;
fnames = fieldnames(TB);
flength = length(fnames);

% Prepare the array
D_0D = ones(3,1000,flength).*nan;
D_2D = ones(3,1000,flength).*nan;
tau_0D = ones(3,1000,flength).*nan;
tau_2D = ones(3,1000,flength).*nan;

for ktest = 1:flength
    tau0D = [];
    tau2D = [];
    D0D   = [];
    D2D   = [];
    if ~isempty(TB.(fnames{ktest}).D)
        t = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).ID.n;
        t0D = TB.(fnames{ktest}).t0D.*TB.(fnames{ktest}).ID.n;
        % Vector stress ==================================================%
        tau0D = TB.(fnames{ktest}).tau0D(3,:); % 0D Numerical experiment
        tau2D = TB.(fnames{ktest}).tau(1:length(t)); % 2D numerical experiment
        % Vector Necking rate
        D0D = TB.(fnames{ktest}).D0D2; % 0D numerical experiment
        D2D   = TB.(fnames{ktest}).D;    % 2D numerical experiment
        %==================================================================%
        D_0D(1,1:length(t),ktest)=t;
        D_0D(2,1:length(t),ktest)=D0D;
        D_0D(3,1:length(t),ktest)= ones(length(t),1)*TB.(fnames{ktest}).res;
        %=================================================================%
        D_2D(1,1:length(t),ktest)=t;
        D_2D(2,1:length(t),ktest)=D2D;
        D_2D(3,1:length(t),ktest)= ones(length(t),1)*TB.(fnames{ktest}).res;
        %=================================================================%
        tau_0D(1,1:length(t0D),ktest)=t0D;
        tau_0D(2,1:length(t0D),ktest)=tau0D;
        tau_0D(3,1:length(t0D),ktest)= ones(length(t0D),1)*TB.(fnames{ktest}).res;
        %=================================================================%
        tau_2D(1,1:length(t),ktest)=t;
        tau_2D(2,1:length(t),ktest)=tau2D;
        tau_2D(3,1:length(t),ktest)= ones(length(t),1)*TB.(fnames{ktest}).res;
    end
end

path2colormap2 = strcat('Utilities\ScientificColourMaps8\','viridis','\','viridis','.mat');
load(path2colormap2);
cmap = colormap(viridis);

fun_0D = Manuscript_function_Container;

c0D = squeeze(D_0D(3,:,:))*100;

color_lists = fun_0D.color_computation(flength,c0D,0,10,0);

f=figure(1);
clf;

set(gcf, 'Units','centimeters', 'Position', [0, 0, 15,15], 'PaperUnits', 'centimeters', 'PaperSize', [15, 15])
subplot(2,2,1)
% Create the axis
ax0 = gca();
p0x  = 0.08;
p0y  = 0.08; 
sz_axh = 0.35;
sz_axv = 0.4;
dh     = sz_axh+0.05;
dv     = sz_axv+0.05; 

ax0.Position = [p0x,p0y+dv,sz_axh,sz_axv];

x = squeeze(D_2D(1,:,:));
y = squeeze(D_2D(2,:,:));
c = squeeze(D_2D(3,:,:));

for i = 1:flength
    hold on
    aa = x(:,i);
    bb = y(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    if ~isempty(aa)
        pl0(i)=plot(ax0,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.9);
    end
    hold off
end
caxis(ax0,[0,10]);
ax0 = set_default_necking_rate(ax0); 
ax0.XTickLabel = [];
ax0.XLabel.String = []; 
ax0.YLabel.Units = 'normalized';
ax0.YLabel.Position(2)=0.50; 
ax0.YLabel.Position(1)=-0.01; 
xline(ax0,1.0,'Color','#501414',LineStyle="-",LineWidth=1.2,Alpha=0.8)



subplot(2,2,2)
ax2 = gca();
ax2.Position = [p0x+dh,p0y+dv,sz_axh,sz_axv];

x = squeeze(D_0D(1,:,:));
y = squeeze(D_0D(2,:,:));
c = squeeze(D_0D(3,:,:));

for i = 1:flength
    hold on
    aa = x(:,i);
    bb = y(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    if ~isempty(aa)
        pl0(i)=plot(ax2,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.9);
    end
    hold off
end
caxis(ax2,[0,10]);
ax2 = set_default_necking_rate(ax2); 
ax2.XTickLabel = [];
ax2.XLabel.String = []; 
ax2.YLabel.String = [];
ax2.YTickLabel = [];
xline(ax2,1.0,'Color','#501414',LineStyle="-",LineWidth=1.2,Alpha=0.8)


% Set default value for necking rate
% Set default
subplot(2,2,3)
ax3 = gca();
ax3.Position = [p0x,p0y,sz_axh,sz_axv];

x = squeeze(tau_0D(1,:,:));
y = squeeze(tau_0D(2,:,:));
c = squeeze(tau_0D(3,:,:));

for i = 1:flength
    hold on
    aa = x(:,i);
    bb = y(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    if ~isempty(aa)
        pl0(i)=plot(ax3,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.9);
    end
    hold off
end
caxis(ax3,[0,10]);
ax3 = set_default_tau(ax3); 
ax3.YLabel.Units = 'normalized';
ax3.YLabel.Position(2) = 0.5;
ax3.YLabel.Position(1) = -0.01;
xline(ax3,1.0,'Color','#501414',LineStyle="-",LineWidth=1.2,Alpha=0.8)
yline(ax3,1.0,'Color','#660000',LineStyle="-",LineWidth=1.2,Alpha=0.8)



subplot(2,2,4)
ax4 = gca();
ax4.Position = [p0x+dh,p0y,sz_axh,sz_axv];

x = squeeze(tau_2D(1,:,:));
y = squeeze(tau_2D(2,:,:));
c = squeeze(tau_2D(3,:,:));

for i = 1:flength
    hold on
    aa = x(:,i);
    bb = y(:,i);
    cc = c(:,i);
    aa = aa(~isnan(aa));
    bb = bb(~isnan(aa));
    cc = cc(~isnan(aa));
    if ~isempty(aa)
        pl3(i)=plot(ax4,aa,bb,'Color',cmap(color_lists(i),:),'LineWidth',0.9);
    end
    hold off
end
caxis(ax4,[0,10]);
ax4 = set_default_tau(ax4); 
ax4.YLabel.String = [];
ax4.YTickLabel = []; 
xline(ax4,1.0,'Color','#501414',LineStyle="-",LineWidth=1.2,Alpha=0.8)
yline(ax4,1.0,'Color','#660000',LineStyle="-",LineWidth=1.2,Alpha=0.8)
cbar = colorbar(ax4,"eastoutside");
cbar.Box = 'on';
cbar.LineWidth = 1.2;
cbar.Color = 'k';
cbar.TickLabelInterpreter = 'latex';
cbar.Label.String = 'Fit,$\%$';
cbar.Label.Interpreter = 'latex';
cbar.Label.Position(1)=0.15;
cbar.Label.Color = 'w';
cbar.Position= [p0x+2*dh,p0y+sz_axv/4,0.05,sz_axv+0.5*sz_axv+0.05];

t0 = annotation("textbox");
t0.Position = [p0x+sz_axh-0.075,p0y+dv+sz_axv-0.05,0.05,0.05];
t0.Interpreter = 'latex';
t0.String = '$[2D]$';
t0.FontSize = 12;
t0.FontWeight = 'bold';

t0.LineStyle = 'none';

t1 = annotation("textbox");
t1.Position = [p0x+dh+sz_axh-0.075,p0y+dv+sz_axv-0.05,0.05,0.05];
t1.Interpreter = 'latex';
t1.String = '$[0D]$';
t1.FontSize = 12;
t1.FontWeight = 'bold';
t1.LineStyle = 'none';
figure_name = 'figure_5.png';
filename = fullfile(ptsave,figure_name);
exportgraphics(f,filename,'Resolution',1000,'BackgroundColor','white')         






end
%
function ax = set_default_necking_rate(ax)

ax.Box = 'on';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1.2;
ax.XLim = [0.1,10];
ax.XTick = [0.1000    1.0000   10.0000];
ax.XTickLabel = {'0.1','1','10'};
ax.XLabel.String = '$t^{\dagger}$';
ax.XTickLabelRotation = 0;
ax.XLabel.Units = 'normalized';
ax.XLabel.Position(2)=-0.1;
ax.XLabel.Interpreter = 'latex';
ax.LineWidth = 1.2;
ax.YColor = [0,0,0];
ax.TickLabelInterpreter = 'latex';
ax.YLabel.String = '$D^{\dagger}$';
ax.YLabel.Interpreter = 'latex';
ax.YTick = [0.1,0.5,1.0];
ax.YLim = [0.1,1.0];
ax.YTickLabel = {'$0.1$',[],'$1$'};
ax.Box ='on';
ax.FontUnits = 'centimeters';
ax.FontSize = 0.5;

end

function ax = set_default_tau(ax)

ax.Box = 'on';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1.2;
ax.XLim = [0.1,10];
ax.XTick = [0.000    1.0000   10.0000];
ax.XTickLabel = {'0','1','10'};
ax.XLabel.String = '$t^{\dagger}$';
ax.XTickLabelRotation = 0;
ax.XLabel.Units = 'normalized';
ax.XLabel.Position(2)=-0.05;
ax.XLabel.Interpreter = 'latex';
ax.LineWidth = 1.2;
ax.YColor = [0,0,0];
ax.TickLabelInterpreter = 'latex';
ax.YLabel.String = '$\tau^{\dagger}$';
ax.YLabel.Interpreter = 'latex';
ax.YTick = [0.1,5.0,10.0];
ax.YLim = [0.0,10.0];
ax.YTickLabel = {'$0.0$',[],'$10.0$'};
ax.Box ='on';
ax.FontUnits = 'centimeters';
ax.FontSize = 0.5;
end