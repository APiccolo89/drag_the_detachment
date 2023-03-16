    %=========================================================================%
    % Realistic data
    %
    %=========================================================================%
    clear all;
    close all;
    addpath Realistic_Case\
    addpath Utilities\
    addpath Adimensional\
    addpath ../../../../freezeColors-master/
    addpath ../../../export_fig/
    load('Data_Base_Realistic_Scenario.mat');
    
    size_marker = 12;
    LineWidth = 1.0;
    font_legend = 8;
    Tp    = 1350+273.15;
    Pr    = 3300*100e3*9.81;
    t0    = 100e6;
    Vnv   = [0:0.1:30].*10^-6;
    Vdv   = Vnv;
    T_mean   = [600:25:1100]+273.15;
    D0 = 80e3;
    L0 = [300,400,500,600].*1e3;
    
    %%
    xiuS = 1e6;
    %% Figure Manuscript Dry Olivine Hirth Kolshdet
    x = squeeze(UPPER_MANTLE.Vd(:,:,1)).*1e6;
    y = squeeze(UPPER_MANTLE.Vn(:,:,1)).*1e6;
    levels_cr = 1:1:10;
    level_la  = 17:1:26;
    vel  = [0,5];
    vel  = vel./100;
    vel  = vel./(365.25.*60.*60.*24);
    ec_m  = vel./1000e3;
    ec   = UM.Bn.*exp(-(UM.En+10e-6.*Pr)./(UM.R.*1100)).*t0.^UM.n;
    
    eps   = ec_m./ec;
    
    figure(1)
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 13, 20], 'PaperUnits', 'centimeters', 'PaperSize', [12, 20])
    
    xiuM = squeeze(UPPER_MANTLE.xiumP(:,:,1));
    eta_eff_REF = (1./squeeze(UPPER_MANTLE.eta0DMP(:,:,1))+1./squeeze(UPPER_MANTLE.eta0MP(:,:,1))).^(-1);
    
    eta_eff_0  = squeeze(UPPER_MANTLE.eta0DMP(:,:,1))./(1+xiuM.*eps(1).^((UM.n-1)./UM.n));
    eta_eff_5  = squeeze(UPPER_MANTLE.eta0DMP(:,:,1))./(1+xiuM.*eps(2).^((UM.n-1)./UM.n));
    
    levels_cr = 1:1:10;
    level_la  = 17:1:26;
    
    tiledlayout(7,4,'TileSpacing','tight', 'Padding', 'none');
    nexttile(1,[2,1]);
    a=pcolor(x,y,log10(xiuM));shading interp;
    colormap(crameri('nuuk',length(levels_cr)-1));
    ylabel('$V_n [10^{-6}\frac{m^3}{J}]$',Interpreter='latex')
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];
    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];
    axis_x.TickLabels = [];
    % title(['$\xi^{UM}, L_0 = ',num2str(L0(i)./1000),' [km]$'],Interpreter='latex')
    caxis([levels_cr(1),levels_cr(end)])
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
    line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
    line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
    %axis square;
    box on
    grid on
    set(gca,'Layer','top')
    scatter(6,15,40,'k','filled','x')
    hold off
    freezeColors;
    nexttile(2,[2,1]);
    c=pcolor(x,y,log10(eta_eff_0));shading interp;
    colormap(crameri('bilbao',length(level_la)-1));
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];
    axis_x.TickLabels   = [];
    
    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];
    set(gca,'yticklabel',[])
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = LineWidth)
    line([10,10],[2,27],'Color', 'r',LineWidth = LineWidth)
    line([2,10],[2,2],'Color' ,'r',LineWidth = LineWidth)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = LineWidth)
    scatter(6,15,size_marker,'k','filled','x')
    hold off
    caxis([level_la(1),level_la(end)])
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    
    nexttile(3,[2,1])
    pcolor(x,y,log10(eta_eff_5));shading interp;
    colormap(crameri('bilbao',length(level_la)-1));
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];
    axis_x.TickLabels   = [];
    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];
    set(gca,'yticklabel',[])
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = LineWidth)
    line([10,10],[2,27],'Color', 'r',LineWidth = LineWidth)
    line([2,10],[2,2],'Color' ,'r',LineWidth = LineWidth)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = LineWidth)
    hold off
    caxis([level_la(1),level_la(end)])
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    
    nexttile(4,[2,1])
    pcolor(x,y,log10(eta_eff_REF));shading interp;
    colormap(crameri('bilbao',length(level_la)-1));
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];
    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];
    set(gca,'yticklabel',[])
    set(gca,'xticklabel',[])
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = LineWidth)
    line([10,10],[2,27],'Color', 'r',LineWidth = LineWidth)
    line([2,10],[2,2],'Color' ,'r',LineWidth = LineWidth)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = LineWidth)
    hold off
    caxis([level_la(1),level_la(end)])
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    
    
    
    
    xiuM = squeeze(UPPER_MANTLE.xiumP(:,:,end));
    eta_eff_REF = (1./squeeze(UPPER_MANTLE.eta0DMP(:,:,end))+1./squeeze(UPPER_MANTLE.eta0MP(:,:,end))).^(-1);
    
    eta_eff_0  = squeeze(UPPER_MANTLE.eta0DMP(:,:,end))./(1+xiuM.*eps(1).^((UM.n-1)./UM.n));
    eta_eff_5  = squeeze(UPPER_MANTLE.eta0DMP(:,:,end))./(1+xiuM.*eps(2).^((UM.n-1)./UM.n));
    
    nexttile(9,[2,1])
    a=pcolor(x,y,log10(xiuM));shading interp;
    colormap(crameri('nuuk',length(levels_cr)-1));
    xlabel('$V_d [10^{-6}\frac{m^3}{J}]$',Interpreter='latex')
    ylabel('$V_n [10^{-6}\frac{m^3}{J}]$',Interpreter='latex')
    
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];
    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];
    caxis([levels_cr(1),levels_cr(end)])
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = LineWidth)
    line([10,10],[2,27],'Color', 'r',LineWidth = LineWidth)
    line([2,10],[2,2],'Color' ,'r',LineWidth = LineWidth)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = LineWidth)
    hold off
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    
    nexttile(10,[2,1])
    
    c=pcolor(x,y,log10(eta_eff_0));shading interp;
    colormap(crameri('bilbao',length(level_la)-1));
    xlabel('$V_d [10^{-6}\frac{m^3}{J}]$',Interpreter='latex')
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];
    
    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];
    
    set(gca,'yticklabel',[])
    
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
    line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
    line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
    hold off
    caxis([level_la(1),level_la(end)])
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    
    nexttile(11,[2,1])
    pcolor(x,y,log10(eta_eff_5));shading interp;
    colormap(crameri('bilbao',length(level_la)-1));
    xlabel('$V_d [10^6\frac{m^3}{J}]$',Interpreter='latex')
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];
    
    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];
    
    set(gca,'yticklabel',[])
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = LineWidth)
    line([10,10],[2,27],'Color', 'r',LineWidth = LineWidth)
    line([2,10],[2,2],'Color' ,'r',LineWidth = LineWidth)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = LineWidth)
    hold off
    caxis([level_la(1),level_la(end)])
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    
    nexttile(12,[2,1])
    pcolor(x,y,log10(eta_eff_REF));shading interp;
    colormap(crameri('bilbao',length(level_la)-1));
    xlabel('$V_d [10^6\frac{m^3}{J}]$',Interpreter='latex')
    
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];
    
    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];
    
    set(gca,'yticklabel',[])
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = LineWidth)
    line([10,10],[2,27],'Color', 'r',LineWidth = LineWidth)
    line([2,10],[2,2],'Color' ,'r',LineWidth = LineWidth)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = LineWidth)
    hold off
    caxis([level_la(1),level_la(end)])
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    nexttile(17,[1,2])
    axis off
    colormap(crameri('nuuk',length(levels_cr)-1));
    caxis([levels_cr(1),levels_cr(end)])
    c=colorbar(gca,"south",TickLabelInterpreter="latex");
    c.Label.Interpreter = 'latex';
    c.Label.String      = '$log_{10}\left(\xi^{UM}\right), [n.d.]$';
    freezeColors;freezeColors(c)
    nexttile(19,[1,2])
    axis off
    axis off
    colormap(crameri('bilbao',length(level_la)-1));
    caxis([level_la(1),level_la(end)])
    c=colorbar(gca,"south",TickLabelInterpreter="latex");
    c.Label.Interpreter = 'latex';
    c.Label.String      = '$log_{10}\left(\eta_{eff}\right), [Pas]$';
    freezeColors;freezeColors(c)
    
    
    
    nexttile(21,[2,4])
    
    % ax_1.Visible = 'off';
    ivd =  find(Vdv.*1e6==6);
    ind =  find(Vdv.*1e6>=2 & Vdv.*1e6<=27);
    ind2 = find(Vdv.*1e6>=2 & Vdv.*1e6<=12);
    ivn =  find(Vnv.*1e6==15);
    a=squeeze(SLAB.MMT_av(ivd,ivn,:)-273.15);
    b=squeeze(SLAB.eta0DS(ivd,ivn,:));
    b_m=squeeze(SLAB.eta0DS(ind2(1),ind(1),:));
    b_M=squeeze(SLAB.eta0DS(ind2(end),ind(end),:));
    
    c=squeeze(SLAB.eta0S(ivd,ivn,:));
    c_m=squeeze(SLAB.eta0S(ind2(1),ind(1),:));
    c_M=squeeze(SLAB.eta0S(ind2(end),ind(end),:));
    eff = (1./b+1./c).^(-1);
    d=squeeze(UPPER_MANTLE.eta0DMP(ivd,ivn,:));
    e=squeeze(UPPER_MANTLE.eta0MP(ivd,ivn,:));
    f=squeeze(UPPER_MANTLE.xiumP(ivd,ivn,:));
    
    hold on
    %p(1)=plot(a,b_m,'LineWidth',LineWidth,'Color','r','LineStyle',':');
    p(2)=plot(a,b,'LineWidth',LineWidth,'LineStyle',':','Color','r');
    %p(3)=plot(a,b_M,'LineWidth',LineWidth,'Color','r','LineStyle',':');
    %p(4)=plot(a,c_m,'LineWidth',LineWidth,'Color','b','LineStyle','--');
    p(5)=plot(a,c,'LineWidth',LineWidth,'LineStyle','--','Color','b');
    %p(6)=plot(a,c_M,'LineWidth',LineWidth,'Color','b','LineStyle','--');
    p(7)=line([600 1100],[1e24 1e24],'Color','k','LineWidth',LineWidth,'LineStyle','-');
    p(8)=line([600 1100],[1e18 1e18],'Color','k','LineWidth',LineWidth,'LineStyle','-');
    patch([a' fliplr(a')],[b_m' fliplr(b_M')],'r','FaceAlpha',0.1);
    patch([a' fliplr(a')],[c_m' fliplr(c_M')],'b','FaceAlpha',0.1);
    
    hold off
    ax_1 = gca;
    
    ax_1.YScale='log';
    ax_1.XColor = [0,0,0];
    ax_1.YColor = [0,0,0];
    ax_1.XLabel.String=['$T ^\circ C$'];
    ax_1.XLabel.Interpreter = 'latex';
    ax_1.YLabel.Interpreter = 'latex';
    ax_1.YLabel.String      = '$\eta_{0} \mathrm{[Pas]}$';
    ax_1.TickLabelInterpreter ='latex';
    ax_1.XGrid ='on';
    ax_1.YGrid ='on';
    ax_1.XMinorGrid ='on';
    ax_1.YMinorGrid ='on';
    ax_1.Box = 'on';
    ax_1.XAxis.LineWidth=LineWidth;
    ax_1.YAxis.LineWidth=LineWidth;
    l= legend('$\eta^S_{0,\mathrm{d}} V_\mathrm{d} = 6[10^{-6}\frac{m^3}{J}]$','$\eta^S_{0,\mathrm{n}},V_\mathrm{n} = 15[10^{-6}\frac{m^3}{J}]$');
    l.Interpreter = 'latex';
    l.FontSize = font_legend;
    l.Location = 'northeast';
    
    F1D.a.Color = 'w';
    filename = 'Realistic_data_Dry_Olivine';
    folder   = '../Manuscript_Main';
    if not(isfolder(folder))
        mkdir(folder);
    end
    pt=fullfile(folder,filename);
    set(gcf, 'Color', 'w')
    export_fig(pt,'-r600')
    %% ======================================================================= %
    % Slab
    %%=========================================================================%
    s0 = 100e6;
    n  = 3.5;
    nlm       = Problem_type.NonLinear;   % Switching the position of linear-non_linear activate the non linear upper mantle routine.
    
    ivd_av =  find(Vdv.*1e6==6);
    ivd_m  =  find(Vdv.*1e6==4);
    ivd_M  =  find(Vdv.*1e6==9);
    ind =  find(Vdv.*1e6>=2 & Vdv.*1e6<=27);
    ind2 = find(Vdv.*1e6>=2 & Vdv.*1e6<=12);
    ivn_av =  find(Vnv.*1e6==15);
    ivn_m =  find(Vnv.*1e6==8);
    ivn_M = find(Vnv.*1e6==20);
    vnv = (2.0:1.0:27.0).*10^-6;
    vdv = (2.0:1.0:12.0).*10^-6;
    ind_n = zeros(length(vnv),1);
    ind_d = zeros(length(vdv),1);
    for i = 1:length(vnv)
        i = i 
        ind_n(i)=find(Vnv==vnv(i));
    end
    for i = 1:length(vdv)
        ind_d(i)=find(Vnv==vnv(i));
    end

  %  ind_n = [ind];
  %  ind_d = [ind2];
    [ind_nm,ind_dm,L0m,T_meanm] = ndgrid(ind_n,ind_d,L0,T_mean);
    xium_      = ind_nm.*0.0;
    eta0D_     = xium_;
    data_n = length(ind_nm(:));
    disp(data_n)
    
    for it = 1:data_n
        T_name   = strcat('T_',num2str(it));
        in = ind_nm(it);
        id = ind_dm(it);
        l0 = L0m(it);
        T_Slab = T_meanm(it);
        eta0DS = SLAB.eta0DS(id,in,T_mean == T_Slab);
        xiUS   = SLAB.xiuS(id,in,T_mean == T_Slab);
        eta0DM = UPPER_MANTLE.eta0DMP(id,in,L0 == l0);
        xiUM   = UPPER_MANTLE.xiumP(id,in,L0 == l0);
    
    
        [Temp]   = Processing_simulation(eta0DS,xiUS,n,l0,t0,D0,eta0DM,0,xiUM,nlm);
    
    
        disp([num2str(it),'out of',num2str(data_n), 'tests'])
        if isreal(Temp.D_norm) && Temp.initial_data.Psi<1
            Tests.(T_name)=Temp;
            Meta_dataReal.Vd = Vdv(id);
            Meta_dataReal.Vn = Vnv(in);
            Meta_dataReal.T_Slab=T_Slab;
            Meta_dataReal.Pr    = Pr;
            Meta_dataReal.Tp    = Tp;
            [Meta_dataReal.Cd,Meta_dataReal.Cn,Meta_dataReal.phi] = UM.Compute_Cd_Cn(Pr,Vnv(in),Vdv(id),Tp);
            Meta_dataReal.w    = UM.rho.*UM.g;
            Tests.(T_name).Meta_dataReal = Meta_dataReal;
            disp(['Detachment_time = ',num2str(Tests.(T_name).t_det*3.5)]);
            disp(['T = ', num2str(T_Slab-273.15)])
            disp(['log10 etaM = ', num2str(log10(eta0DM))])
            disp(['log10 etaS = ', num2str(log10(eta0DS))])
            disp(['log10 xiUM = ', num2str(log10(xiUM))])
            disp(['log10 xiUS = ', num2str(log10(xiUS))])
            disp(['Vd         = ', num2str(Vdv(id)*1e6)])
            disp(['Vn = ', num2str((Vnv(in)*1e6))])
            disp(['l0 = ', num2str((l0./1e3))])
        end
    end
%%
[Data_S] = extract_information_detachment(Tests,0,nlm);
save('Realistic_DataSet.mat','Tests','Data_S')




