
function plot_the_figures(S,pt_save,nlm)


FIT = S.FIT;

% Figure 1 [Lambda fit, and fetch 1 parameter]

Figure1 = scatter_plot_post_process;
Figure1.figure_number = 1;
Figure1.x = FIT.Lambda;
Figure1.y = FIT.fetch(1,:);
Figure1.c = FIT.xium;
Figure1.logx = 'log';
Figure1.logy = 'linear';
Figure1.logcolor = 1;
Figure1.colormap_f ='oslo';
Figure1.xlabel = '$\Lambda_0$';
Figure1.ylabel = '$\omega_{depth}$';
Figure1.clabel = [];
Figure1.xlim = [min(FIT.Lambda)-0.01*mean(FIT.Lambda),1.0];
Figure1.ylim = [0,0.45];

Figure1.save_path = pt_save;
Figure1.name_figure = 'Lambda_fetch1';
Figure1.size_picture = [12.5,13];
Figure1.make_scatter_plot;

Figure2 = scatter_plot_post_process;
Figure2.figure_number = 2;
Figure2.x = FIT.Lambda;
Figure2.y = FIT.fetch(2,:);
Figure2.c = FIT.xium;
Figure2.logx = 'log';
Figure2.logy = 'linear';
Figure2.logcolor = 1;
Figure2.colormap_f ='oslo';
Figure2.xlabel = '$\frac{\Lambda}{1+\xi^{M}\tau^{M,(n-1)}}$';
Figure2.ylabel = '$\omega_{s}$';
Figure2.xlim = [min(FIT.Lambda)-0.01*mean(FIT.Lambda),1.0];
Figure2.clabel = [];
Figure2.xlim = [min(FIT.Lambda)-0.01*mean(FIT.Lambda),1.0];
Figure2.save_path = pt_save;
Figure2.name_figure = 'Lambda_fetch2';
Figure2.size_picture = [12.5,13];
Figure2.make_scatter_plot;


%%

Figure3 = scatter_plot_post_process;
Figure3.figure_number = 3;
Figure3.x = FIT.Lambda;
Figure3.y = FIT.fitting_p.*100;
Figure3.c = [102,102,255]./255;
Figure3.logx = 'log';
Figure3.logy = 'linear';
Figure3.logcolor = 1;
Figure3.colormap_f =[];
Figure3.xlabel = '$\frac{\Lambda}{1+\xi^{M}\tau^{M,(n-1)}}$';
Figure3.xlim = [min(FIT.Lambda)-0.01*mean(FIT.Lambda),1.0];
Figure3.ylabel = '$F [\%]$';
Figure3.clabel = [];
Figure3.xlim = [min(FIT.Lambda)-0.01*mean(FIT.Lambda),1.0];
Figure3.ylim = [0,28];

Figure3.save_path = pt_save;
Figure3.name_figure = 'Lambda_fitting';
Figure3.size_picture = [12.5,13];
Figure3.make_scatter_plot;

%%


TB = S.TB;
% 1D line plot
fnames = fieldnames(TB);
flength = length(fnames);

%Figure Stress
%Figure Thickness
%Figure Topography
%Figure Topography n.d.
F_S = ones(3,1000,flength).*nan;
F_Sft = ones(3,1000,flength).*nan;
F_D = ones(3,1000,flength).*nan;
F_Dft = ones(3,1000,flength).*nan;
F_T = ones(3,1000,flength).*nan;
F_Tnd = ones(3,1000,flength).*nan;
F_dTdt = ones(3,1000,flength).*nan;
F_dTdt_nd =ones(3,1000,flength).*nan;




for ktest = 1:flength
    if ~isempty(TB.(fnames{ktest}).D)
        t = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).ID.n;
        if nlm.Linear == 0
            Lambda = (TB.(fnames{ktest}).ID.Lambda)./(1+TB.(fnames{ktest}).ID.Df_UM.*TB.(fnames{ktest}).ID.tau_mc.^(3.5-1));
        else
            Lambda = (TB.(fnames{ktest}).ID.Lambda);
        end


        % Fill the vector with dimensionless arrays: [2D numerical test]
        %=====================================================================%
        F_S(1,1:length(t),ktest)= t;
        F_D(1,1:length(t),ktest)= t;
        F_Tnd(1,1:length(t),ktest)= t;
        F_dTdt_nd(1,1:length(t),ktest)= t;
        %=====================================================================%
        % Fill the vector with dimensional time
        t_dim = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).tc;
        t_dim = t_dim./(365.25*60*60*24*1e6);
        %=====================================================================%
        F_T(1,1:length(t),ktest)= t_dim;
        F_dTdt(1,1:length(t),ktest)= t_dim;
        %=====================================================================%
        % Fill the fitting time
        t0D = TB.(fnames{ktest}).t0D.*TB.(fnames{ktest}).ID.n;
        %======================================================================%
        F_Sft(1,1:length(t0D),ktest)= t0D;
        F_Dft(1,1:length(t),ktest)= t;
        %Fill the properties
        D = TB.(fnames{ktest}).D;
        D_0D = TB.(fnames{ktest}).D0D2;
        tau = TB.(fnames{ktest}).tau;
        tau0D = TB.(fnames{ktest}).tau0D;
        %=====================================================================%
        F_S(2,1:length(t),ktest) = tau(1:length(t));
        F_Sft(2,1:length(t0D),ktest) = tau0D(3,:);
        F_D(2,1:length(t),ktest) = D(1:length(t));
        F_Dft(2,1:length(t),ktest) = D_0D;
        %=====================================================================%
        %=====================================================================%
        %Fill the color
        F_S(3,1:length(t),ktest) =ones(length(t),1).*TB.(fnames{ktest}).res.*100;
        F_Sft(3,1:length(t0D),ktest) =ones(length(t0D),1).*TB.(fnames{ktest}).res.*100;
        F_D(3,1:length(t),ktest) =ones(length(t),1).*TB.(fnames{ktest}).res.*100;
        F_Dft(3,1:length(t),ktest) =ones(length(t),1).*TB.(fnames{ktest}).res.*100;

        %=====================================================================%

        topo = TB.(fnames{ktest}).topo;
        dt_dim = t_dim(2:1:end)-t_dim(1:1:end-1);
        dTopo = topo(2:1:length(t_dim))-topo(1:1:length(t_dim)-1);
        dU = dTopo./dt_dim;
        T_Mean = (t_dim(2:1:end)+t_dim(1:end-1)).*0.5;
        t_M = 0.5.*(t(1:1:end-1)+t(2:1:end));
        %=====================================================================%
        F_Tnd(1,1:length(t),ktest)= t;
        F_dTdt_nd(1,1:length(t_M),ktest)= t_M;
        F_Tnd(2,1:length(t),ktest)= topo(1:length(t))./(-min(topo(1:length(t))));
        F_dTdt_nd(2,1:length(t_M),ktest)= dU;
        %
        F_T(1,1:length(t_dim),ktest)= t_dim;
        F_dTdt(1,1:length(T_Mean),ktest)= T_Mean;
        F_T(2,1:length(t_dim),ktest)= topo(1:length(t_dim))./(-min(topo(1:length(t))));
        F_dTdt(2,1:length(T_Mean),ktest)= dU;
        %
        F_T(3,1:length(T_Mean),ktest)= log10(ones(length(T_Mean),1).*Lambda);
        F_dTdt(3,1:length(T_Mean),ktest)= log10(ones(length(T_Mean),1).*Lambda);
        F_Tnd(3,1:length(T_Mean),ktest)= log10(ones(length(T_Mean),1).*Lambda);
        F_dTdt_nd(3,1:length(T_Mean),ktest)= log10(ones(length(T_Mean),1).*Lambda);
    end
end
path2colormap = strcat('Utilities\ScientificColourMaps8\','lipari','\','lipari','.mat');
load(path2colormap);
cmap2 = colormap(lipari);

path2colormap2 = strcat('Utilities\ScientificColourMaps8\','glasgow','\','glasgow','.mat');
load(path2colormap2);
cmap = colormap(glasgow);

fun_0D = Manuscript_function_Container;  

c0D = squeeze(F_Dft(3,:,:));
c2D = squeeze(F_D(3,:,:));

color_lists_FIT0D = fun_0D.color_computation(flength,c0D,0,20,0);
color_lists_FIT2D = fun_0D.color_computation(flength,c2D,0,20,0);



Figure4 = line_plot_post_process;
Figure4.figure_number = 4;
Figure4.logx = 'linear';
Figure4.logy = 'linear';
Figure4.colormap_f = cmap;
Figure4.xlabel = '$t^{\dagger}$';
Figure4.ylabel = '$D^{\dagger,2D}$';
Figure4.clabel = '$F[\%]$';
Figure4.clim  = [0,20];
Figure4.ctick = [0:5:20];
Figure4.xlim  = [0.01,10];
Figure4.ylim  = [0.1,1];
Figure4.x      = squeeze(F_D(1,:,:));
Figure4.y      = squeeze(F_D(2,:,:));
Figure4.c      = log10(squeeze(F_D(3,:,:)));
Figure4.size_picture = [13,12.5];
Figure4.name_figure = 'Thickness_2D';
Figure4.save_path = pt_save;
Figure4.multiline_plot(flength,color_lists_FIT2D)

Figure5 = line_plot_post_process;
Figure5.figure_number = 5;
Figure5.logx = 'linear';
Figure5.logy = 'linear';
Figure5.colormap_f = cmap;
Figure5.xlabel = '$t^{\dagger}$';
Figure5.ylabel = '$D^{\dagger,0D}$';
Figure5.clabel = '$F[\%]$';
Figure5.clim  = [0,20];
Figure5.ctick = [0:5:20];
Figure5.xlim  = [0.01,10];
Figure5.ylim  = [0.1,1.0];
Figure5.x      = squeeze(F_Dft(1,:,:));
Figure5.y      = squeeze(F_Dft(2,:,:));
Figure5.c      = log10(squeeze(F_Dft(3,:,:)));
Figure5.size_picture = [13,12.5];
Figure5.name_figure = 'Thickness_0D';
Figure5.save_path = pt_save;
Figure5.multiline_plot(flength,color_lists_FIT0D)


cL = 10.^squeeze(F_T(3,:,:));
color_lists_Lambda = fun_0D.color_computation(flength,cL,-4,0);



Figure6 = line_plot_post_process;
Figure6.figure_number = 6;
Figure6.logx = 'linear';
Figure6.logy = 'linear';
Figure6.colormap_f = cmap2;
Figure6.xlabel = '$t^{\dagger}$';
Figure6.ylabel = '$H [km]$';
Figure6.clabel = '$\frac{\Lambda_c}{1+\xi^M\tau_0^{M,n-1}}$';
Figure6.clim  = [-4,0];
Figure6.ctick = [-4:1:0];
Figure6.xlim  = [0.01,10];
Figure6.ylim  = [-1.1,0.0];
Figure6.x      = squeeze(F_Tnd(1,:,:));
Figure6.y      = squeeze(F_Tnd(2,:,:));
Figure6.c      = log10(squeeze(F_Tnd(3,:,:)));
Figure6.size_picture = [13,12.5];
Figure6.name_figure = '2D_2fetch_topography_d';
Figure6.save_path = pt_save;
Figure6.multiline_plot(flength,color_lists_Lambda)

Figure7 = line_plot_post_process;
Figure7.figure_number = 7;
Figure7.logx = 'linear';
Figure7.logy = 'linear';
Figure7.colormap_f = cmap2;
Figure7.xlabel = '$t [Myrs]$';
Figure7.ylabel = '$H [km]$';
Figure7.clabel = '$\frac{\Lambda_c}{1+\xi^M\tau_0^{M,n-1}}$';
Figure7.clim  = [-4,0];
Figure7.ctick = [-4:1:0];
Figure7.xlim  = [0.0,35.0];
Figure7.ylim  = [-1.1,0.0];
Figure7.x      = squeeze(F_T(1,:,:));
Figure7.y      = squeeze(F_T(2,:,:));
Figure7.c      = log10(squeeze(F_T(3,:,:)));
Figure7.size_picture = [13,12.5];
Figure7.name_figure = '2D_2fetch_topography_dim_d';
Figure7.save_path = pt_save;
Figure7.multiline_plot(flength,color_lists_Lambda)

Figure8 = line_plot_post_process;
Figure8.figure_number = 8;
Figure8.logx = 'linear';
Figure8.logy = 'linear';
Figure8.colormap_f = cmap2;
Figure8.xlabel = '$t^{\dagger}$';
Figure8.ylabel = '$\dot{H} [km/yrs]$';
Figure8.clabel = '$\frac{\Lambda_c}{1+\xi^M\tau_0^{M,n-1}}$';
Figure8.clim  = [-4,0];
Figure8.ctick = [-4:1:0];
Figure8.xlim  = [0.01,10];
Figure8.ylim   = [-0.05,0.2];
Figure8.x      = squeeze(F_dTdt_nd(1,:,:));
Figure8.y      = squeeze(F_dTdt_nd(2,:,:));
Figure8.c      = log10(squeeze(F_dTdt_nd(3,:,:)));
Figure8.size_picture = [13,12.5];
Figure8.name_figure = '2D_2fetch_topography';
Figure8.save_path = pt_save;
Figure8.multiline_plot(flength,color_lists_Lambda)


Figure9 = line_plot_post_process;
Figure9.figure_number = 9;
Figure9.logx = 'linear';
Figure9.logy = 'linear';
Figure9.colormap_f = cmap2;
Figure9.xlabel = '$t [Myrs]$';
Figure9.ylabel = '$\dot{H} [km/yrs]$';
Figure9.clabel = '$\frac{\Lambda_c}{1+\xi^M\tau_0^{M,n-1}}$';
Figure9.clim  = [-4,0];
Figure9.ctick = [-4:1:0];
Figure9.x      = squeeze(F_dTdt(1,:,:));
Figure9.y      = squeeze(F_dTdt(2,:,:));
Figure9.c      = log10(squeeze(F_dTdt(3,:,:)));
Figure9.xlim   = [0.1,35];
Figure9.ylim   = [-0.05,0.2];
Figure9.size_picture = [13,12.5];
Figure9.name_figure = '2D_2fetch_topography_dim';
Figure9.save_path = pt_save;
Figure9.multiline_plot(flength,color_lists_Lambda)


end
