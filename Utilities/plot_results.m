
function plot_results(Tests,name,ptsave,nlm)
    % Input 
    % Testdata=> Data Structure containing all the tests 
    % Short description
    % Plot t/tc - Dnorm of all the tests and colored as a function of
    % lambda.
    % Plot td against Lambda scatter plot.
    % Stress => To Do How to retrieve the stress data per each timestep 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Data_S] = extract_information_detachment(Tests,0,nlm);
    % function plot not dimensional data {+ additional collect data for the scatter plot
    % i.e. Lambda, t_d,Psi} 
    Plot_1D_Plots(Tests,Data_S,name,ptsave,nlm);
    % function to do scatter plot
    %plot_scatter(Data_S,name,ptsave,'t_det',nlm)
    %plot_scatter(Data_S,name,ptsave,'time_max-time det',nlm)
    
end

function [Data_S] = Plot_1D_Plots(Tests,Data_S,name,ptsave,pb)

    tic 
    %plot3D_setExp(Tests,Data_S,name,'D_norm',ptsave,pb)
    toc 
    tic
    %plot1D_setExp(Tests,Data_S,name,'D_norm',ptsave,pb)
    toc
    tic
    %plot1D_setExp(Tests,Data_S,name,'tau_eff',ptsave,pb)
    toc
    tic
    %plot1D_setExp(Tests,Data_S,name,'tau_D_tau_B',ptsave,pb)
    toc
    tic
    %plot1D_setExp(Tests,Data_S,name,'epsilon',ptsave,pb)
    toc
    tic
    plot1D_setExp(Tests,Data_S,name,'dDdt_tauD',ptsave,pb)
    toc
    tic
%    plot1D_setExp(Tests,Data_S,name,'dDdt_tauD_norm',ptsave,pb)
    toc
    tic
    plot1D_setExp(Tests,Data_S,name,'dDdt_t_d',ptsave,pb)
    toc
    if pb.islinear == 0
        plot1D_setExp(Tests,Data_S,name,'Lambda',ptsave,pb)
        plot1D_setExp(Tests,Data_S,name,'eta_mantle',ptsave,pb)

        tic
    end
end

function plot_scatter(Data_S,name,ptsave,field,nlm) 
c = Data_S(9,:);
double = 0; 
figure(1)

if strcmp(field,'t_det')
    x = Data_S(1,:);
    if nlm.islinear ==0
        x = x./(Data_S(13,:));
    end
    y = 1.0./Data_S(3,:);
    ylabel_ = ('$\frac{n}{t^{det}_c}$ [n.d]');
    ylim_=[0,0.0];
    fin = 'Global_test1DS_time_det';
elseif strcmp(field,'tau_ii')
    x = Data_S(1,:);
    y = Data_S(4,:); % tau max
    y2 = Data_S(6,:); % tau @ detachment
    ylabel_=('$\frac{\tau_{eff}}{\tau_{B_0}}$');
    ylim_ = [1,7.0];
    double = 1.0 ; 
    fin = 'Global_test1DS_Stress';

elseif strcmp(field,'time_max-time det')
    x = Data_S(1,:);
    y = Data_S(3,:)-Data_S(5,:);
    ylabel_=('$t_{det}-t(\tau_{eff}^{MAX}) [n.d.]$');
    fin = 'Global_test1DS_dT';
    ylim_ = [0,0];
elseif strcmp(field,'t_det_error')
    x = Data_S(1,:);
    y = Data_S(7,:);
    ylabel_ =('t^D_{det}-t^A_{det} [n.d.]');
    fin = 'Global_test1_error';
    ylim_ = [0,0];


elseif strcmp(field,'error_gl')
    x = Data_S(1,:);
    y = Data_S(8,:);
    y_err = Data_S(9,:)-Data_S(10,:);
    ylabel_ = ('$mean(D_A(t)-D_D(t)$, mean error dimension D, adimensional D ');
    fin = 'Global_test1_error_GL';
    ylim_ = [0,0];

end
s = 1;%Data_S(9,:);

%scatter(x(Data_S(8,:)==0.5),y(Data_S(8,:)==0.5),10,"black",'filled','d')

hold on 
scatter(x(Data_S(8,:)==100.0),y(Data_S(8,:)==100.0),10,c,'filled','o')


%if double >0 
 %   scatter(x,y2,5,c,'filled','o','MarkerEdgeColor','k','LineWidth',0.5)
  %  legend({'$max(\frac{\tau_{II}}{\tau_{B,0}})$','$\frac{\tau^{Det}_{II}}{\tau_{B,0}}$'},'Interpreter','latex');
 %end

try
    cmap = colormap(crameri('nuuk'));
catch
    cmap = colormap('jet');
end
% Set colorbar
c=colorbar;
c.Label.String = 'log10(\Psi) [n.d.]';
grid on
box on
xlabel('$log_{10}(\frac{\Lambda}{\xi^{UM}})$[n.d.]',Interpreter='latex')
ylabel(ylabel_,Interpreter='latex')
if ylim_(2) == 0
    disp('none')
else
    ylim(ylim_);
end
xlim([10^(-8),1])
%ylabel('t^O_d/t^P_d [n.d]')
set(gca, 'XScale', 'log')
xlim([10^(-7.3),10^(0)])
name_picture = strcat(fin,name,'.png');
pt=fullfile(ptsave,'SCAT');
if not(isfolder(pt))
     mkdir(pt);
end
title(name)
pt=fullfile(pt,name_picture);

print(pt,'-dpng')
clf; 
close;

end 
function  plot1D_setExp(Tests, Data_S,name,field,ptsave,nlm)


    % Pre Process Data 
    %=====================================================================%
    % Interpolate in a  common ground all the data (i.e. vector between  %
    % 0-100
    %======================================================================
    z_min = -6;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0; %round(log10(max(Data_S(1,Data_S(1,:)<1.0)))); %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
    figure(1)
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    
    if strcmp(field,'D_norm')
        VAL = 0;
    elseif strcmp(field,'tau_eff')
        VAL=1.0; 
    elseif strcmp(field,'epsilon')
        VAL = 2.0;
    elseif strcmp(field,'tau_D_tau_B')
        VAL = 3.0;
    elseif strcmp(field,'Lambda')
        VAL =4.0;
    elseif strcmp(field,'eta_mantle')
        VAL = 5.0; 
    elseif strcmp (field,'dDdt_tauD')
        VAL = 6.0; 
    elseif strcmp(field,'dDdt_t_d')
        VAL = 7.0; 
    else 
        VAL = 5; 
    end

    try
        cmap = colormap(crameri('bilbao'));
    catch
        cmap = colormap('jet');
    end
    % Set colorbar
    % Set colorbar
   
    i = 1; % iterator 
    fn = fieldnames(Tests);

    for k = 1:numel(fn) 
       TD = Tests.(fn{k});
       if ~isnan(Data_S.tdet(k))
           if VAL == 0 
               buf = TD.D_norm; 
               ylim([0.1,1.0])
               ylabel('$\frac{D}{D_0} [n.d.]$','interpreter','latex')
           elseif VAL == 1
               buf = TD.tau(3,:);
               ylabel('$\frac{\tau_{eff}}{\tau_{B,0}} [n.d.]$','interpreter','latex')
           elseif VAL==2 
               buf = TD.eps(1,:);
               ylabel('$\frac{\dot{\varepsilon}_{II}}{\dot{\varepsilon}_{B,0}} [n.d.]$','interpreter','latex')
               set(gca, 'YScale', 'log')
           elseif VAL ==3 || VAL == 5
               buf = (-TD.tau(2,:)./TD.tau(1,:));

               ylabel('$\frac{\tau_{D}}{\tau_{B}} [n.d.]$','interpreter','latex')
               
               set(gca, 'YScale', 'log')
           elseif VAL == 4
               buf = (TD.Lambda);
               ylabel('$ \Lambda(t) [n.d.]$','interpreter','latex')
                set(gca, 'YScale', 'log')
           elseif VAL==5.0
               buf = (TD.Lambda./TD.initial_data.Lambda).*TD.initial_data.eta0DM;
               ylabel('$ \eta^{UM}_{eff}(t) [n.d.]$','interpreter','latex')
               set(gca, 'YScale', 'log')

           elseif VAL==6 
               buf = (-TD.tau(2,:)./TD.tau(1,:));
               buf = (buf(1:1:end-1)+buf(2:1:end))/2; 
               buf = buf./max(abs(buf));
           else 
                
    
           end
      
           hold on 
           if VAL  <5
               x = TD.time*TD.initial_data.n; 
               xlabel('$\frac{t}{t_d} [n.d.]$','interpreter','latex')
               if strcmp(field,'D_norm')
                x = (TD.time*TD.initial_data.tc)./(365.25*60*60*24*1e6); 
                xlabel('${t} [Myrs]$','interpreter','latex')
               else
                   xlim([0,15])
               end

           else
               xlabel('$log_{10}\left(\frac{dD}{dt}\right) [n.d.]$','Interpreter','latex')
               ylabel('$log_{10}\left(\frac{\tau_{D}}{\tau_{0,B}}\right) [n.d.]$','Interpreter','latex')

               dD = TD.D_norm; 
               dD = dD(2:1:end)-dD(1:1:end-1);
               x = TD.time*TD.initial_data.n;  
               x2 = x(2:1:end)/2+x(1:1:end-1)/2;
               time_max=TD.time_t_M*TD.initial_data.n;
               ind = find(x2<=time_max);
               ind = ind(end)-1; 
               x = x(2:1:end)-x(1:1:end-1);
               x = dD./x; 
               if VAL == 7
                buf = x; 
                x   = x2; 
                ylabel('$log_{10}\left(\frac{dD}{dt}\right) [n.d.]$','Interpreter','latex')
                xlabel('$log_{10}\left(\frac{t}{t_d}\right) [n.d.]$','Interpreter','latex')


               end
               %if VAL == 7 
                %    xlabel('$norm\left(\frac{dD}{dt})\right) [n.d.]$','Interpreter','latex')
                 %   ylabel('$norm\left(\frac{\tau_D}{\tau_{0,B}}\right) [n.d.]$','Interpreter','latex')
                  %  x = x./max(abs(x)); 
               
               %end
               scatter(x(ind),buf(ind),10,'red','filled','d');
               
           end
           

           plot(x,buf,'Color','k',LineWidth=0.7)
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')
       end        
       i = i+1; 
    end
   
    box on
    box on
    grid on

    name_picture = strcat('Global_test1D',field,name,'.png');
    pt=fullfile(ptsave,field);
    if not(isfolder(pt))
         mkdir(pt);
    end
    title(name)

    pt=fullfile(pt,name_picture);

    print(pt,'-dpng')
    clf; 
    close;
end

function  plot3D_setExp(Tests, Data_S,name,field,ptsave,nlm)

close all; 
clf; 
% Pre Process Data
%=====================================================================%
% Interpolate in a  common ground all the data (i.e. vector between  %
% 0-100
%======================================================================
if nlm.islinear
    z_min = -6;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0; %round(log10(max(Data_S(1,Data_S(1,:)<1.0)))); %log10(max(Data_S(1,:)));
else
    z_min = -12;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0; 
end
% See if the user installed Crameri cmap utilities, otherwise punish
% him with jet colormap by default
% Shamelessly copied from
% https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis

fn = fieldnames(Tests);

t_intp = 0:0.0001:20;
DS = ones(3000,length(fn)).*nan;
ts = DS; 
bla = 0; 
for k = 1:numel(fn)
    TD = Tests.(fn{k});
    y_d = TD.D_norm;
    is = isreal(y_d);
    if is == 0 
        bla = bla+1 ; 
    end
    if is == 1 
        %y_d = real(y_d); 
        x = TD.time.*TD.initial_data.n;
        DS(1:length(y_d),k) = y_d;
        ts(1:length(y_d),k) = x; 
        if nlm.islinear ==0
            lambda=log10(TD.initial_data.Lambda./(1+TD.initial_data.Df_UM));
            V = (log10(TD.initial_data.Lambda./(1+TD.initial_data.Df_UM))-z_min)/(z_max-z_min);
        else
            lambda = log10(TD.initial_data.Lambda);
            V = (log10(TD.initial_data.Lambda)-z_min)/(z_max-z_min);
        end
        if V<0
            V=0;
        elseif V>1
            V=1;
        end
        V=round(1+V*(256-1));%round to nearest index
        V_map(k) = V;
        L_axis(k) = 10.^lambda;
    
        Det(k) = x(end);
    else
        V_map(k) = nan;
        L_axis(k)=nan;
        DS(:,k) = nan;
        ts(:,k) = nan;
        Det(k) = nan;    
    end
end


disp(bla)


F=figure(1);
% Collect the field names
% Set the min and max of lambda value for the coloring of the plot
% Set colorbar
% Set colorbar
 try
        cmap = colormap(crameri('bilbao'));
    catch
        cmap = colormap('jet');
    end
c=colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = ['$log_{10}\left(\frac{\Lambda}{1+\xi^{UM}}\right), [n.d.]$'];
% Define coloraxis
caxis([z_min z_max])
% min/max are rounded using the log values of Lambda
Ticks=round((z_min)):round((z_max));
% Produce the tick
TickLabels=arrayfun(@(x) sprintf('10^{%d}',x),Ticks,'UniformOutput',false);
try
    set(c,'Ticks',Ticks);
    set(c,'TickLabels',TickLabels);
catch %HG1
    TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
    set(c,'YTick',Ticks);
    set(c,'YTickLabel',TickLabels);
end

V_map(isnan(V_map))=1; 
C = cmap(V_map,:);
bla = 1; 
for i = 1:(size(DS,2))
    hold on
    plot3(ts(:,i),DS(:,i),L_axis(i).*ones(length(ts(:,1)),1),'Color',C(i,:),'LineWidth',1.2);
    
end
scatter3(Det(~isnan(Det)==1),0.05.*ones(length(Det(~isnan(Det)==1)),1),L_axis(~isnan(Det)==1),10,'k','filled','d')
plot3(ts(:,1:end),DS(:,1:end),1e-14.*ones(length(ts(:,1)),1),'Color','k','LineWidth',0.6)
set(gca, 'ZScale', 'log');
set(gca,'color',[0.7 0.7 0.7])
xlim([0,15])
ylim([0.05,1.1])
grid on 
box on
view(-40,30)
xlabel('$\frac{t}{t_d}, [n.d]$',Interpreter='latex')
ylabel('$\frac{D}{D_0}, [n.d]$',Interpreter='latex')
zlabel('$log_{10}\left(\frac{\Lambda}{1+\xi^{UM}}\right), [n.d.]$',Interpreter='latex')


name_picture = strcat('3dDnorm.png');
pt=fullfile(ptsave,field);
if not(isfolder(pt))
    mkdir(pt);
end
pt=fullfile(pt,name_picture);
disp(pt)
set(F, 'InvertHardCopy', 'off');
set(F, 'PaperPositionMode','auto', 'PaperOrientation','portrait')
export_fig(F,pt, '-dpng','-r600','-transparent')
clf;
close;
end


function  plot3D_setExp2(Tests, Data_S,name,field,ptsave,nlm)


% Pre Process Data
%=====================================================================%
% Interpolate in a  common ground all the data (i.e. vector between  %
% 0-100
%======================================================================
z_min = -6;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
z_max = 0; %round(log10(max(Data_S(1,Data_S(1,:)<1.0)))); %log10(max(Data_S(1,:)));
% See if the user installed Crameri cmap utilities, otherwise punish
% him with jet colormap by default
% Shamelessly copied from
% https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis

fn = fieldnames(Tests);

t_intp = 0:0.0001:20;
DS = ones(3000,length(fn)).*nan;
ts = DS; 
bla = 0; 
for k = 1:numel(fn)
    TD = Tests.(fn{k});
    y_d = (-TD.tau(2,:)./TD.tau(1,:));
    is = isreal(y_d);
    if is == 0 
        bla = bla+1 ; 
    end
    if is == 1 
        %y_d = real(y_d); 
        x = TD.time.*TD.initial_data.n;
        DS(1:length(y_d),k) = y_d;
        ts(1:length(y_d),k) = x; 
        if nlm.islinear ==0
            lambda=log10(TD.initial_data.Lambda./(1+TD.initial_data.Df_UM));
            V = (log10(TD.initial_data.Lambda./(1+TD.initial_data.Df_UM))-z_min)/(z_max-z_min);
        else
            lambda = log10(TD.initial_data.Lambda);
            V = (log10(TD.initial_data.Lambda)-z_min)/(z_max-z_min);
        end
        if V<0
            V=0;
        elseif V>1
            V=1;
        end
        V=round(1+V*(256-1));%round to nearest index
        V_map(k) = V;
        L_axis(k) = 10.^lambda;
    
        Det(k) = x(end);
    else
        V_map(k) = nan;
        L_axis(k)=nan;
        DS(:,k) = nan;
        ts(:,k) = nan;
        Det(k) = nan;    
    end
end


disp(bla)


F=figure(1);
% Collect the field names
% Set the min and max of lambda value for the coloring of the plot
% Set colorbar
% Set colorbar
 try
        cmap = colormap(crameri('bilbao'));
    catch
        cmap = colormap('jet');
    end
c=colorbar;
c.Label.String = 'log10(\Lambda) [n.d.]';
% Define coloraxis
caxis([z_min z_max])
% min/max are rounded using the log values of Lambda
Ticks=round((z_min)):round((z_max));
% Produce the tick
TickLabels=arrayfun(@(x) sprintf('10^{%d}',x),Ticks,'UniformOutput',false);
try
    set(c,'Ticks',Ticks);
    set(c,'TickLabels',TickLabels);
catch %HG1
    TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
    set(c,'YTick',Ticks);
    set(c,'YTickLabel',TickLabels);
end

V_map(isnan(V_map))=1; 
C = cmap(V_map,:);
bla = 1; 
for i = 1:(size(DS,2))
    hold on
    plot3(ts(:,i),DS(:,i),L_axis(i).*ones(length(ts(:,1)),1),'Color',C(i,:),'LineWidth',1.2);
    
end
plot3(ts(:,1:end),DS(:,1:end),1e-10.*ones(length(ts(:,1)),1),'Color','k','LineWidth',0.6)
set(gca, 'ZScale', 'log');
set(gca, 'YScale', 'log');

set(gca,'color',[0.7 0.7 0.7])
xlim([0,15])
ylim([0.05,1.1])
grid on 
box on
view(-40,30)
xlabel('$\frac{t}{t_d}, [n.d]$',Interpreter='latex')
ylabel('$\frac{D}{D_0}, [n.d]$',Interpreter='latex')
zlabel('$log_{10}(\Lambda), [n.d.]$',Interpreter='latex')


name_picture = strcat('3D_other.png');
pt=fullfile(ptsave,field);
if not(isfolder(pt))
    mkdir(pt);
end
pt=fullfile(pt,name_picture);
disp(pt)
set(F, 'InvertHardCopy', 'off');
set(F, 'PaperPositionMode','auto', 'PaperOrientation','portrait')
export_fig(F,pt, '-dpng','-r600','-transparent')
clf;
close;
end