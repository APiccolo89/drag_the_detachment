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
nlm = Problem_type.NonLinear;
Df_UM = nan;
DB_path = '../Test_DB2.hdf5';
pt_save = '../Results/Comparison2D_0DNL';
if not(isfolder(pt_save))
     mkdir(pt_save);
end
l = h5info(DB_path,'/Viscous/HR_DIS');
%l2 =h5info(DB_path,'/Viscous/LR_NLM');
TestsA  = {l.Groups.Name};
%TestsB  = {l2.Groups.Name};
HR = 1:length(TestsA);
%LR = -(1:length(TestsB));
Res = [HR]; 
Tests   = [TestsA]; 
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
    ylabel('$\omega_{\ell}$',Interpreter='latex')
    grid on
    box on
    set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')

    pt=fullfile(pt_save,'FS3');
    print(pt,'-dpng')
 figure(9)
    title_plot = ('DET');
    scatter((FIT.Lambda(FIT.Res>0)),1./FIT.Detachment(1,FIT.Res>0),"k",'filled','o')
       hold on
    scatter((FIT.Lambda(FIT.Res<0)),1./FIT.Detachment(1,FIT.Res<0),"k",'+')
   
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$\omega_{\ell}$',Interpreter='latex')
    grid on
    box on
    set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')
    ylim([0,1])
    xlim([10^-8,10^-0])
    pt=fullfile(pt_save,'FS4');
    print(pt,'-dpng')

     figure(10)
    title_plot = ('DET_REAL');
    scatter((FIT.Lambda(FIT.Res>0)),1./FIT.Detachment(2,FIT.Res>0),"k",'filled','o')  
      hold on
    scatter((FIT.Lambda(FIT.Res<0)),1./FIT.Detachment(2,FIT.Res<0),"k",'+')
   
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$\omega_{\ell}$',Interpreter='latex')
    grid on
    box on
    ylim([0,1])
    xlim([10^-6,10^-1])
    set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')

    pt=fullfile(pt_save,'FS5');
    print(pt,'-dpng')


