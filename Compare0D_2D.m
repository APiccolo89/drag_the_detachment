%=========================================================================%

clear all;
close all;
addpath Adimensional\
addpath Dimensional\
addpath Utilities\
addpath D0_D2_comparison_function\
clear all
% Common Data
D0=80e3;
n=3.5;
Df_S=10;
nlm = Problem_type.Linear;
Df_UM = nan;
DB_path = '../Test_DB2.hdf5';
pt_save = '../Results/Comparison2D_0D';
if not(isfolder(pt_save))
     mkdir(pt_save);
end
l = h5info(DB_path,'/Viscous/HR');
l2 =h5info(DB_path,'/Viscous/LR_NLM');
TestsA  = {l.Groups.Name};
TestsB  = {l2.Groups.Name};
HR = 1:length(TestsA);
LR = -(1:length(TestsB));
Res = [HR LR]; 
Tests   = [TestsA TestsB]; 
%function Optimize Data Base 

[TB,FIT] = perform_optimization_DataBase(Tests,n,D0,Df_S,nlm,Df_UM,DB_path,pt_save,Res);

clf; 
close all; 

    figure(6)
    title_plot = ('$log_{10}(\Lambda)$ vs $fit$');
    scatter((FIT.Lambda(FIT.Res>0)),FIT.fitting_p(FIT.Res>0).*100,"k",'filled','o')
    hold on
    scatter((FIT.Lambda(FIT.Res<0)),FIT.fitting_p(FIT.Res<0).*100,"k",'+')
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$F \%$',Interpreter='latex')
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale','log')
    grid on
    box on
    pt=fullfile(pt_save,'FS1');
    print(pt,'-dpng')

    figure(7)
    title_plot = ('Effective dDdt');
    scatter((FIT.Lambda(FIT.Res>0)),FIT.fetch(2,FIT.Res>0),"k",'filled','o')
    hold on
    scatter((FIT.Lambda(FIT.Res<0)),FIT.fetch(2,FIT.Res<0),"k",'+')
   
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$\omega_{\frac{dD}{dt}}$',Interpreter='latex')
    grid on
    box on
    set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')

    pt=fullfile(pt_save,'FS2');
    print(pt,'-dpng')

    figure(8)
    title_plot = ('Effective dDdt');
    scatter((FIT.Lambda(FIT.Res>0)),FIT.fetch(1,FIT.Res>0),"k",'filled','o')
    hold on
    scatter((FIT.Lambda(FIT.Res<0)),FIT.fetch(1,FIT.Res<0),"k",'+')
   
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$\omega_{\frac{dD}{dt}}$',Interpreter='latex')
    grid on
    box on
    set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')

    pt=fullfile(pt_save,'FS3');
    print(pt,'-dpng')


    figure(9)
    title_plot = ('Effective dDdt');
    scatter((FIT.Lambda(FIT.Res>0)),abs(FIT.fetch(2,FIT.Res>0)),"k",'filled','o')
    hold on
    scatter((FIT.Lambda(FIT.Res<0)),abs(FIT.fetch(2,FIT.Res<0)),"k",'+')
   
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$\omega_{\frac{dD}{dt}}$',Interpreter='latex')
    grid on
    box on
    set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')

    pt=fullfile(pt_save,'FS41');
    print(pt,'-dpng')


  figure(9)
    title_plot = ('Effective dDdt');
    scatter(FIT.fetch(1,FIT.Res>0),FIT.fetch(2,FIT.Res>0),"k",'filled','o')
    hold on
    scatter(FIT.fetch(1,FIT.Res<0),FIT.fetch(2,FIT.Res<0),"k",'+')
   
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$\omega_{\frac{dD}{dt}}$',Interpreter='latex')
    grid on
    box on
    set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')

    pt=fullfile(pt_save,'FS4');
    print(pt,'-dpng')

    x = (FIT.Lambda(FIT.Res>0));
    x2= (FIT.Lambda(FIT.Res<0));
    A = 1-abs(FIT.fetch(1,FIT.Res>0));
    B = 1-abs(FIT.Dp(FIT.Res>0)./2);
    y = (B-A)./B; 

    A2 = abs(FIT.fetch(1,FIT.Res<0));
    B2 = abs(FIT.Dp(FIT.Res<0)./2);
    y2 = (B2-A2)./B2;  
    figure(10)
    %title_plot = ('$log_{10}(\Lambda)$ vs $n_{eff}$');
    scatter(x,abs(y).*100,"k",'filled','o')
    hold on 
    %scatter((Lambda),0.5-((T(:,1)-T(:,3)))./2,'r','filled')
    scatter(x2,abs(y2).*100,"k",'+')
    ylim([0,50])
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$\omega_l$',Interpreter='latex')
    grid on
    box on
    set(gca, 'XScale', 'log')
    pt=fullfile(pt_save,'FS4');
    print(pt,'-dpng')


    figure(11)
    %title_plot = ('$log_{10}(\Lambda)$ vs $n_{eff}$');
    scatter((FIT.Lambda),FIT.Detachment(1,:),"k",'filled','o')
    hold on 
    %scatter((FIT.Lambda),FIT.Detachment(2,:),'b','filled','square')
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$\frac{t^{2D}_{det}-t^{0D}_{det}}{\frac{t^{2D}_{det}}$',Interpreter='latex')
    grid on
    box on 
    set(gca, 'XScale', 'log')
    pt=fullfile(pt_save,'FS5');
    print(pt,'-dpng')


