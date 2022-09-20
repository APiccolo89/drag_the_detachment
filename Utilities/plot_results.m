
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
       Data_S(3,i) = TD.initial_data.n*TD.t_det;
       Data_S(4,i) = TD.t_t_max;
       Data_S(5,i) = TD.time_t_M*TD.initial_data.n;
       Data_S(6,i) = TD.t_t_det;
       Data_S(7,i) = abs(TD.t_det-TD.Testdata_a.t_det); 
       Data_S(8,i) = TD.Interpolation.error_DA(1);
       Data_S(9,i) = TD.Interpolation.error_DA(2);
       Data_S(10,i) = TD.Interpolation.error_DA(3);
       i = i+1;
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
end

function plot_scatter(Data_S,name,ptsave,field) 
c = Data_S(2,:);
double = 0; 
if strcmp(field,'t_det')
    x = Data_S(1,:);
    y = Data_S(3,:);
    ylabel('t^O_d/t^P_d [n.d]')
    fin = 'Global_test1DS_time_det';
elseif strcmp(field,'tau_ii')
    x = Data_S(1,:);
    y = Data_S(4,:); % tau max
    y2 = Data_S(6,:); % tau @ detachment
    ylabel('$\frac{\tau_{eff}}{\tau_{B_0}}$','Interpreter','latex');
    double = 1.0 ; 
    fin = 'Global_test1DS_Stress';

elseif strcmp(field,'time_max-time det')
    x = Data_S(1,:);
    y = Data_S(3,:)-Data_S(5,:);
    ylabel('$t_{det}-t(\tau_{eff}^{MAX}) [n.d.]$','Interpreter','latex');
    fin = 'Global_test1DS_dT';

elseif strcmp(field,'t_det_error')
    x = Data_S(1,:);
    y = Data_S(7,:);
    ylabel('t^D_{det}-t^A_{det} [n.d.]');
    fin = 'Global_test1_error';


elseif strcmp(field,'error_gl')
    x = Data_S(1,:);
    y = Data_S(8,:);
    y_err = Data_S(9,:)-Data_S(10,:);
    ylabel('mean(D_A(t)-D_D(t), mean error dimension D, adimensional D ','Interpreter','latex');
    fin = 'Global_test1_error_GL';
end
figure(1)
scatter(x,y,20,c,'filled','d')
hold on 

if double >0 
    scatter(x,y2,5,c,'filled','o')
%elseif strcmp(field,'error_gl')
 %   errorbar(x,y,y_err,'|r');
end

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
xlabel('\Lambda [n.d.]')
%ylabel('t^O_d/t^P_d [n.d]')
%ylim([0.8,20])
if strcmp(field,'t_det')
    ylabel('t^O_d/t^P_d [n.d]')
 elseif strcmp(field,'tau_ii')
    ylabel('$\frac{\tau_{eff}}{\tau_{B_0}}$','Interpreter','latex'); 
elseif strcmp(field,'time_max-time det')
    ylabel('$t_{det}-t(\tau_{eff}^{MAX}) [n.d.]$','Interpreter','latex');
 end
set(gca, 'XScale', 'log')
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
    z_min = log10(min(Data_S(1,:)));
    z_max  =log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
    if strcmp(field,'D_norm')
        VAL = 0;
    elseif strcmp(field,'tau_eff')
        VAL=1; 
    elseif strcmp(field,'epsilon')
        VAL = 2.0
    else 
        VAL = 3.0; 
    end

    try
        cmap = colormap(crameri('Bilbao'));
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
       else
           buf = log10(-TD.tau(2,:)./TD.tau(1,:));
           ylabel('$\frac{\tau_{D}}{\tau_{B}} [n.d.]$','interpreter','latex')
           set(gca, 'YScale', 'log')
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
       plot(TD.time*TD.initial_data.n,buf,'Color',C)
       grid on 
       xlim([0,10])
       xlabel('$n*(\frac{t}{t_c}) [n.d.]$','interpreter','latex')
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