
function plot_results(Tests,name,ptsave)
    % Input 
    % Testdata=> Data Structure containing all the tests 
    % Short description
    % Plot t/tc - Dnorm of all the tests and colored as a function of
    % lambda.
    % Plot td against Lambda scatter plot.
    % Stress => To Do How to retrieve the stress data per each timestep 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % function plot not dimensional data {+ additional collect data for the scatter plot
    % i.e. Lambda, t_d,Psi} 
    [Data_S] =  Plot_1D_Plots(Tests,name,ptsave);
    % function to do scatter plot
    plot_scatter(Data_S,name,ptsave,'t_det');
    plot_scatter(Data_S,name,ptsave,'t_det_error');
    plot_scatter(Data_S,name,ptsave,'error_gl');
    plot_scatter(Data_S,name,ptsave,'tau_ii');
    plot_scatter(Data_S,name,ptsave,'time_max-time det');



end

function [Data_S] = Plot_1D_Plots(Tests,name,ptsave)

    % Collect the field names 
    fn = fieldnames(Tests);
    % number of test 
    itest = length(fn);
    % Prepare Data_S array
    Data_S = zeros(10,itest); 
    i = 1;
    for k = 1:numel(fn) 
     
       TD = Tests.(fn{k});
       
       Data_S(1,i) = TD.initial_data.Lambda;
       Data_S(2,i) = log10(TD.initial_data.Psi); 
       if ~isempty(TD.t_det) && isreal(TD.D_norm)
        Data_S(3,i) = TD.initial_data.n*TD.t_det;
        Data_S(4,i) = TD.t_t_max;
        Data_S(5,i) = TD.time_t_M*TD.initial_data.n;
        Data_S(6,i) = TD.t_t_det;
        Data_S(7,i) = abs(TD.t_det-TD.Testdata_a.t_det);
        Data_S(10,i) = 1.0 ; 
       else
        Data_S(3,i) = nan;
        Data_S(4,i) = nan;
        Data_S(5,i) = nan;
        Data_S(6,i) = nan;
        Data_S(7,i) = nan;
        Data_S(10,i) = -1.0 ; 
       end 
       Data_S(8,i) = TD.initial_data.Df_S;
       Data_S(9,i) = TD.initial_data.D0/TD.initial_data.l0;
       i = i+1
    end

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
    plot1D_setExp(Tests,Data_S,name,'Lambda',ptsave)
    tic
    plot1D_setExp(Tests,Data_S,name,'Exp',ptsave)
    toc

end

function plot_scatter(Data_S,name,ptsave,field) 
c = Data_S(2,:);
double = 0; 
figure(1)

if strcmp(field,'t_det')
    x = Data_S(1,:);
    y = Data_S(3,:);
    ylabel_ = ('$\frac{t}{t_d}$ [n.d]');
    ylim_=[1,10];
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
    ylabel_=('$t_{det}-t(\tau_{eff}^{MAX}) [n.d.]$')
    fin = 'Global_test1DS_dT';
    ylim_ = [0,0]
elseif strcmp(field,'t_det_error')
    x = Data_S(1,:);
    y = Data_S(7,:);
    ylabel_ =('t^D_{det}-t^A_{det} [n.d.]');
    fin = 'Global_test1_error';
    ylim_ = [0,0]


elseif strcmp(field,'error_gl')
    x = Data_S(1,:);
    y = Data_S(8,:);
    y_err = Data_S(9,:)-Data_S(10,:);
    ylabel_ = ('$mean(D_A(t)-D_D(t)$, mean error dimension D, adimensional D ');
    fin = 'Global_test1_error_GL';
    ylim_ = [0,0];

end
s = 1;%Data_S(9,:);

scatter(x(Data_S(8,:)==0.5),y(Data_S(8,:)==0.5),10,"black",'filled','d')

hold on 
scatter(x(Data_S(8,:)==1.0),y(Data_S(8,:)==1.0),20,"black",'filled','o')
scatter(x(Data_S(8,:)==10.0),y(Data_S(8,:)==10.0),30,"black",'+')
legend({'$ \xi = 0.1 $','$ \xi = 1.0$','$ \xi = 10.0 $'},'Interpreter','latex');

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
xlabel('$log_{10}(\Lambda)$[n.d.]',Interpreter='latex')
ylabel(ylabel_,Interpreter='latex')
if ylim_(2) == 0
    disp('none')
else
    ylim(ylim_);
end
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
    z_min =  round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = round(log10(max(Data_S(1,Data_S(1,:)<1.0)))); %log10(max(Data_S(1,:)));
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
        cmap = colormap(crameri('hawaii'));
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
       if Data_S(10,i) > 0.0
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
           plot(x,buf,'Color',C,LineWidth=1.3)
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