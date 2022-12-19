
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
    plot1D_setExp(Tests,Data_S,name,'D_norm',ptsave)
    toc
    tic
    plot1D_setExp(Tests,Data_S,name,'tau_eff',ptsave)
    toc
    tic
    plot1D_setExp(Tests,Data_S,name,'tau_D_tau_B',ptsave)
    toc
    tic
    plot1D_setExp(Tests,Data_S,name,'epsilon',ptsave)
    toc
    if pb.islinear == 0
        plot1D_setExp(Tests,Data_S,name,'Lambda',ptsave)
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
function  plot1D_setExp(Tests, Data_S,name,field,ptsave)
    figure(1)
    % Collect the field names 
    fn = fieldnames(Tests);
    % Set the min and max of lambda value for the coloring of the plot
    z_min = -8;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0; %round(log10(max(Data_S(1,Data_S(1,:)<1.0)))); %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
    if strcmp(field,'D_norm')
        VAL = 0;
    elseif strcmp(field,'tau_eff')
        VAL=1; 
    elseif strcmp(field,'epsilon')
        VAL = 2.0;
    elseif strcmp(field,'tau_D_tau_B')
        VAL = 3.0;
    elseif strcmp(field,'Lambda')
        VAL =4.0;
    else 
        VAL = 5; 
    end

    try
        cmap = colormap(crameri('broc'));
    catch
        cmap = colormap('jet');
    end
    % Set colorbar
    % Set colorbar
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
    i = 1; % iterator 
    for k = 1:numel(fn) 
       TD = Tests.(fn{k});
       if ~isnan(Data_S.tdet)
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
           else 
                
    
           end
    
           % Normalize Lambda value w.r.t. the limit that I assumed to be
           % likely
           V = (log10(TD.initial_data.Lambda)-z_min)/(z_max-z_min);
           if V<0 
               V=0;
           elseif V>1
                V=1;
           end
           V=round(1+V*(size(cmap,1)-1));%round to nearest index
           C = cmap(V,:);
           hold on 
           if VAL  <5
               x = TD.time*TD.initial_data.n; 
               xlim([0,10])
               %set(gca, 'XScale', 'log')

               xlabel('$n*(\frac{t}{t_c}) [n.d.]$','interpreter','latex')
           else
               x = TD.D_norm;
               xlim([0.1,1.0])
               xlabel('$\frac{D}{D_0} [n.d.]$','Interpreter','latex')
           end
           plot(x,buf,'Color',C,LineWidth=0.8)
           grid on 
       end        
       i = i+1; 
    end
    box on
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