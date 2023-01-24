function [TB,FIT] = perform_optimization_DataBase(Tests,n,D0,Df_S,nlm,Df_UM,DB_path,pt_save,Res)

suc = 1;
for ktest=1:1:length(Tests)
    time_A = cputime;
    Chosen = Tests{ktest};
    [P_Var,D,t,Tau,d,tdet,topo,tauL] = Reading_Data_Base(Chosen,DB_path,Df_S,nlm);
    [ID] = Compute_slab_characteristics(P_Var.eta0DS,Df_S,n,P_Var.L0,P_Var.s0,D0,P_Var.eta0DM,P_Var.xiUM,nlm);
    % Example
    if P_Var.failed == 0
        [TestData_W] = Run_Simulation_DragA(ID.ID_A,nlm);
        ID.ID_A.cut_off_Mantle = 1.0;
        ID.ID_A.cut_off_Slab   = 1.0;
        %ID.ID_A.fetch(1) = abs(double(-d./2));
        f = [abs(double(-d./2)),1.0];
        fun      = @(x) optimizing_fit(ID,nlm,x,D,t);
        f    = fminsearch(fun,f);
        ID.ID_A.fetch(1) = abs(f(1));
        ID.ID_A.fetch(2) = abs(f(2));
        [TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
        [D,D_0D2D] = clean_data(TestData,D,t);
        [D,D_0D2DW] = clean_data(TestData_W,D,t);

        res = goodnessOfFit(D,D_0D2D,'NMSE');
        disp(['The test @ node ', Chosen, ' show a ',num2str(res*100), ' % fit'])
        
        if nlm.islinear == 0
            FIT.Lambda(suc)=ID.ID_A.Lambda./P_Var.xiUM;
        else
          FIT.Lambda(suc)=ID.ID_A.Lambda;
        end
        FIT.fitting_p(suc)=res;
        FIT.fetch(1,suc)=f(1);
        FIT.fetch(2,suc)=f(2);
        FIT.eta0DM(suc)=P_Var.eta0DM;
        FIT.L0(suc)= P_Var.L0;
        FIT.s0(suc)=P_Var.s0;
        FIT.Dp(suc) = d; 
        FIT.Detachment(1,suc) = n*TestData.t_det;
        FIT.Detachment(2,suc) = tdet; 
        FIT.Res(suc)= Res(ktest);
        % Figure Data Base 
    if Res(ktest)>0 
        ls = '-';
    else
        ls = ':';
    end

    FIG1=figure(1);
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = -8;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0;  %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
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
    V = ( log10(FIT.Lambda(suc))-z_min)/(z_max-z_min);
    if V<0 
          V=0;
    elseif V>1
          V=1;
    end
    V=round(1+V*(size(cmap,1)-1));%round to nearest index
    C = cmap(V,:);
    hold on 
    plot(3.5*TestData.time,TestData.D_norm,'Color',C,LineWidth=0.6,LineStyle=ls);
    grid on 
    box on
    xlim([0,10])
    ylim([0.1,1.05])
    pt=fullfile(pt_save);
    if not(isfolder(pt))
         mkdir(pt);
    end
   




    FIG2=figure(2);
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = -8;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0;  %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
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
    V = ( log10(FIT.Lambda(suc))-z_min)/(z_max-z_min);
    if V<0 
          V=0;
    elseif V>1
          V=1;
    end
    V=round(1+V*(size(cmap,1)-1));%round to nearest index
    C = cmap(V,:);
    hold on 
    plot(3.5.*TestData.time,TestData.tau(3,:),'Color',C,LineWidth=0.6,LineStyle=ls)
    grid on 
    box on
    xlim([0,10])
    ylim([0.1,10.01])
    pt=fullfile(pt_save);
    if not(isfolder(pt))
         mkdir(pt);
    end
    
    
    
    FIG3=figure(3);
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = -8;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0;  %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
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
    V = ( log10(FIT.Lambda(suc))-z_min)/(z_max-z_min);
    if V<0 
          V=0;
    elseif V>1
          V=1;
    end
    V=round(1+V*(size(cmap,1)-1));%round to nearest index
    C = cmap(V,:);
    hold on 

    error = abs((((D_0D2D)./(D))))-1;
    error2 = abs((((D_0D2DW)./(D))))-1;

    p=plot(3.5*t./(tdet),error,'Color',C,LineWidth=0.4,LineStyle=ls);
    alpha(p,.1)
    grid on 
    box on
    xlim([0,10])
    %set(gca,'YScale','log')
    xlim([0.0,1.1])
    ylim([-0.4,0.8])

    pt=fullfile(pt_save);
    if not(isfolder(pt))
         mkdir(pt);
    end

FIG4=figure(4);
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = -8;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0;  %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
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
    V = ( log10(FIT.Lambda(suc))-z_min)/(z_max-z_min);
    if V<0 
          V=0;
    elseif V>1
          V=1;
    end
    V=round(1+V*(size(cmap,1)-1));%round to nearest index
    C = cmap(V,:);
    hold on 
    plot(3.5.*t,D./(D(1)),'Color',C,LineWidth=0.6,LineStyle=ls)
    grid on 
    box on
    xlim([0,10])
    ylim([0.1,1.05])
    pt=fullfile(pt_save);
    if not(isfolder(pt))
         mkdir(pt);
    end
   

FIG5=figure(5);
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = -8;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0;  %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
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
    V = ( log10(FIT.Lambda(suc))-z_min)/(z_max-z_min);
    if V<0 
          V=0;
    elseif V>1
          V=1;
    end
    V=round(1+V*(size(cmap,1)-1));%round to nearest index
    C = cmap(V,:);
    hold on 
    plot(3.5.*t,topo./(-min(topo)),'Color',C,LineWidth=0.6,LineStyle=ls)
    grid on 
    box on
    xlim([0,10])
    ylim([-1.05,0.1])
    pt=fullfile(pt_save);
    if not(isfolder(pt))
         mkdir(pt);
    end
   

FIG6=figure(6);
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = -8;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0;  %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
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
    V = ( log10(FIT.Lambda(suc))-z_min)/(z_max-z_min);
    if V<0 
          V=0;
    elseif V>1
          V=1;
    end
    V=round(1+V*(size(cmap,1)-1));%round to nearest index
    C = cmap(V,:);
    hold on 
    plot(3.5.*t,tauL,'Color',C,LineWidth=0.6,LineStyle=ls)
    grid on 
    box on
    xlim([0,10])
    ylim([-0.01,1.05])
    pt=fullfile(pt_save);
    if not(isfolder(pt))
         mkdir(pt);
    end
   
    FIG7=figure(7);
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = -8;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = 0;  %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
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
    V = ( log10(FIT.Lambda(suc))-z_min)/(z_max-z_min);
    if V<0 
          V=0;
    elseif V>1
          V=1;
    end
    V=round(1+V*(size(cmap,1)-1));%round to nearest index
    C = cmap(V,:);
    hold on 
    plot(3.5.*t,Tau,'Color',C,LineWidth=0.6,LineStyle=ls)
    grid on 
    box on
    xlim([0,10])
    ylim([.5,10.0])
    pt=fullfile(pt_save);
    if not(isfolder(pt))
         mkdir(pt);
    end
   





    t = [];
    D_0D2D = [];
    D_0D2DW = [];
    Tau02D = []; 
    error = [];
    suc = suc+1; 

    else
        disp('=====================================================================')
        disp(['You cannot optimize a failure, but, i can give you motivational tips:'])
        disp([ ' it is not always good being consistent in failing'])
        disp('So long and thanks for all the fish')
        disp('=====================================================================')
        f = nan;
        t = [];
        D = [];
        D_0D2D = [];
    end
    TB.(strcat('T',num2str(ktest))).ID = ID.ID_A;
    TB.(strcat('T',num2str(ktest))).P_Var = P_Var;
    TB.(strcat('T',num2str(ktest))).D = D;
    TB.(strcat('T',num2str(ktest))).t = t;
    TB.(strcat('T',num2str(ktest))).D0D2 = D_0D2D;
    TB.(strcat('T',num2str(ktest))).f = f;
    TB.(strcat('T',num2str(ktest))).res = res;

    B = cputime; 
    time_ktest = (B-time_A)/60;

    disp(['The test @ node ', Chosen, ' took ', num2str(time_ktest),' minutes'])

end
pt=fullfile(pt_save,'Global_T');
print(FIG1,pt,'-dpng')
pt=fullfile(pt_save,'Global_Tau');
print(FIG2,pt,'-dpng')
pt=fullfile(pt_save,'Global_Er');
print(FIG3,pt,'-dpng')
pt=fullfile(pt_save,'Real_Tests_D');
print(FIG4,pt,'-dpng')
pt=fullfile(pt_save,'Real_Tests_tauL');
print(FIG5,pt,'-dpng')
pt=fullfile(pt_save,'Real_Tests_Topo');
print(FIG6,pt,'-dpng')
pt=fullfile(pt_save,'Real_Tests_Tau');
print(FIG7,pt,'-dpng')

%close all; 
%clf; 
end