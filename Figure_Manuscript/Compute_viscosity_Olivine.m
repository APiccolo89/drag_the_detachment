%=========================================================================%
% Depth dependent viscosity calculator.
% ========================================================================%
% [INPUT]: -> [t0, age, Tp, D0,L0,(Vn,Vd)]=> Most difficult parameter     |
% 1st Construct a function that compute the reference viscosities for a   |
% given temperature of reference, and pressure of reference.              |
% 2nd Plot the reference viscosities corrected with the pressure          |
% [OUTPUT]:-> [eta0DM,eta0M,eta0DS,eta0S] = function                      |
%=========================================================================%
clear all
close all
addpath ../Realistic_Case/
close all
set(0,'defaultTextInterpreter','latex'); %trying to set the default




% Reference values:
Tp    = 1400+273.15;
Pr    = 3300*100e3*9.81;
t0    = nan;
Vnv   = [0:1:30].*1e-6;
Vdv   = [0:1:30].*1e-6;
T_mean   = [600:25:1100]+273.15;
%T_mean = [700:50:1100]+273.15; 
age = NaN; 
D0 = 80e3;
L0 = [300,400,500,600].*1e3;
S_initial_data = struct('Tp',Tp,'Pr',Pr,'s0',t0,'Vdv',Vdv,'Vnv',Vnv,'T_mean',T_mean,'D0',D0,'L0',L0); 

if ~isnan(t0)
    name_data_base = (['Data_Base_Realistic_Scenario_',num2str(round(Tp-273.15)),'_',num2str(int(t0)),'.mat']);
else
    name_data_base = (['Data_Base_Realistic_Scenario_',num2str(round(Tp-273.15)),'_stress_L0_3D.mat']);
end

folder = '../Data_Base/Realistic_Data_Base/';
figure_folder = 'figure_4__ALT';
supplementary_folder = 'figure_4__ALT'; 
folder_save=fullfile(folder,figure_folder);
folder_supplementary=fullfile(folder,supplementary_folder);
filename = fullfile(folder,name_data_base);

if not(isdir(folder_save))
    mkdir(folder_save);
end

if not(isdir(folder_supplementary))
    mkdir(folder_supplementary);
end


B_n = 1.1e5;
B_n = correct_data(1.1e5,3.5,1,0,1.0,10e3);
B_d = correct_data(1.5e9,1.0,1.0,3,1.0,10e3);
%B_d = B_d.*10^(-6.0*3.5);
B_d_w =  correct_data(1.0e6,1.0,1.0,3,1000.0,10e3);
B_n_w = correct_data(1600,3.5,1.2,0,1000.0,1);
UM        = Mantle_Unit_Properties(3300,3e-5,1050,B_d,B_n,375e3,530e3,3.5);
UM2       = Mantle_Unit_Properties(3300,3e-5,1050,B_d_w,B_n_w,335e3,520e3,3.5);
S         = Mantle_Unit_Properties(3360,3e-5,1050,B_d,B_n,375e3,530e3,3.5);

% Dry Olivine Data:
[UPPER_MANTLE,SLAB] =  main_function_Real(t0, T_mean, Tp,Pr, D0,L0,Vnv,Vdv,UM,S,age);
% Wet Olivine Data:
[UPPER_MANTLE2,~] =  main_function_Real(t0, T_mean, Tp,Pr, D0,L0,Vnv,Vdv,UM,S,age);

%% 
% Dry and Wet database 
save((filename), 'UPPER_MANTLE2','SLAB','UPPER_MANTLE','UM2','UM','S_initial_data');



font_axes = 16; 
font_legend = 14; 
font_text   = 5; 
size_shit = [13,13.5];
LineWidth = 1.0; 
marker_size = 6;

%%

% for i = 1:length(L0)
%     clf; close all;
%     figure(i)
%     set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%     xiuM = squeeze(UPPER_MANTLE.xiumP(:,:,i));
%     eta_eff_REF = (1./squeeze(UPPER_MANTLE.eta0DMP(:,:,i))+1./squeeze(UPPER_MANTLE.eta0MP(:,:,i))).^(-1);
% 
%     eta_eff_0  = squeeze(UPPER_MANTLE.eta0DMP(:,:,i))./(1+xiuM.*eps(1).^((UM2.n-1)./UM2.n));
%     eta_eff_5  = squeeze(UPPER_MANTLE.eta0DMP(:,:,i))./(1+xiuM.*eps(2).^((UM2.n-1)./UM2.n));
% 
%     a=1./(1./SLAB.eta0DS+1./SLAB.eta0S);
%     ind1 = find(T_mean==600+273.15,1); 
%     ind2 = find(T_mean==800+273.15,1);
%     ind3 = find(T_mean==1000+273.15,1);
%     Psi1 =  (eta_eff_REF)./(squeeze(a(:,:,ind1)));
%     Psi2 =  (eta_eff_REF)./(squeeze(a(:,:,ind2)));
%     Psi3 =  (eta_eff_REF)./(squeeze(a(:,:,ind3)));
% 
% 
% 
%     levels_cr = -6:1:26;
%     level_la  = 14:1:26;
%     level_Psi = -20:1:0; 
% 
%     
%     a=pcolor(Vdv.*1e6,Vnv.*1e6,log10(xiuM));shading interp;
%     colormap(crameri('roma',length(levels_cr)-1));
%     axis_x= get(gca, 'XAxis');
%     axis_x.TickLabelInterpreter = 'latex';
%     axis_x.TickValues   = [5:5:25];
% 
%     axis_y= get(gca, 'YAxis');
%     axis_y.TickLabelInterpreter = 'latex';
%     axis_y.TickValues   = [5:5:25];
% 
% 
%     title(['$\xi^{UM}, L_0 = ',num2str(L0(i)./1000),' [km]$'],Interpreter='latex')
%     caxis([levels_cr(1),levels_cr(end)])
%     hold on
%     line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
%     line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
%     line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
%     line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
%     hold off
%     axis square;
%     box on
%     grid on
%     %set(gca,'Layer','top')
%         filename = (['xium',num2str(i)]);
% 
%     pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%       pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%     print(pt,'-dpng')
% 
%     figure(i)
%     set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%     c=pcolor(Vdv.*1e6,Vnv.*1e6,log10(eta_eff_REF));shading interp;
%     colormap(crameri('bilbao',length(level_la)-1));
%     xlabel('$V_d [10^6\frac{m^3}{J}]$',Interpreter='latex')
%     ylabel('$V_n [10^6\frac{m^3}{J}]$',Interpreter='latex')
%      axis_x= get(gca, 'XAxis');
%     axis_x.TickLabelInterpreter = 'latex';
%     axis_x.TickValues   = [5:5:25];
% 
%     axis_y= get(gca, 'YAxis');
%     axis_y.TickLabelInterpreter = 'latex';
%     axis_y.TickValues   = [5:5:25];
% 
% 
%     hold on
%     line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
%     line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
%     line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
%     line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
%     hold off
%     caxis([level_la(1),level_la(end)])
%     axis square;
%     box on
%     grid on
%     set(gca,'Layer','top')
%      filename = (['eta',num2str(i)]);
%          pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%     print(pt,'-dpng')
% 
% 
%     figure(i)
%     set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%     c=pcolor(Vdv.*1e6,Vnv.*1e6,log10(Psi1));shading interp;
%     title(['$T = $',num2str(T_mean(ind1)-273.15), '$[^{\circ}C]$'],Interpreter='latex')
% 
%     colormap(crameri('berlin',length(level_Psi)-1));
%  xlabel('$V_d [10^6\frac{m^3}{J}]$',Interpreter='latex')
%     ylabel('$V_n [10^6\frac{m^3}{J}]$',Interpreter='latex')
%      axis_x= get(gca, 'XAxis');
%     axis_x.TickLabelInterpreter = 'latex';
%     axis_x.TickValues   = [5:5:25];
% 
%     axis_y= get(gca, 'YAxis');
%     axis_y.TickLabelInterpreter = 'latex';
%     axis_y.TickValues   = [5:5:25];
% 
% 
%     hold on
%     line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
%     line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
%     line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
%     line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
%     hold off
%     caxis([level_Psi(1),level_Psi(end)])
%     axis square;
%     box on
%     grid on
%     set(gca,'Layer','top')
%      filename = (['Psi',num2str(i)]);
%          pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%     print(pt,'-dpng')
%     
%     figure(i)
%     set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%     c=pcolor(Vdv.*1e6,Vnv.*1e6,log10(Psi2));shading interp;
%     colormap(crameri('berlin',length(level_la)-1));
%         title(['$T = $',num2str(T_mean(ind2)-273.15), '$[^{\circ}C]$'],Interpreter='latex')
% 
%  xlabel('$V_d [10^6\frac{m^3}{J}]$',Interpreter='latex')
%     ylabel('$V_n [10^6\frac{m^3}{J}]$',Interpreter='latex')
%      axis_x= get(gca, 'XAxis');
%     axis_x.TickLabelInterpreter = 'latex';
%     axis_x.TickValues   = [5:5:25];
%     colormap(crameri('berlin',length(level_Psi)-1));
% 
%     axis_y= get(gca, 'YAxis');
%     axis_y.TickLabelInterpreter = 'latex';
%     axis_y.TickValues   = [5:5:25];
% 
% 
%     hold on
%     line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
%     line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
%     line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
%     line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
%     hold off
%     caxis([level_Psi(1),level_Psi(end)])
%     axis square;
%     box on
%     grid on
%     set(gca,'Layer','top')
%      filename = (['Psi2',num2str(i)]);
%          pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%     print(pt,'-dpng')
% 
% 
%     figure(i)
%     set(gcf, 'Units','centimeters', 'Position', [0, 0, size_shit(1),size_shit(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_shit(1), size_shit(2)])
% 
%     c=pcolor(Vdv.*1e6,Vnv.*1e6,log10(Psi3));shading interp;
%     colormap(crameri('berlin',length(level_Psi)-1));
%     title(['$T = $',num2str(T_mean(ind3)-273.15), '$[^{\circ}C]$'],Interpreter='latex')
% 
%  xlabel('$V_d [10^6\frac{m^3}{J}]$',Interpreter='latex')
%     ylabel('$V_n [10^6\frac{m^3}{J}]$',Interpreter='latex')
%      axis_x= get(gca, 'XAxis');
%     axis_x.TickLabelInterpreter = 'latex';
%     axis_x.TickValues   = [5:5:25];
% 
%     axis_y= get(gca, 'YAxis');
%     axis_y.TickLabelInterpreter = 'latex';
%     axis_y.TickValues   = [5:5:25];
% 
% 
%     hold on
%     line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
%     line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
%     line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
%     line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
%     hold off
%     caxis([level_Psi(1),level_Psi(end)])
%     axis square;
%     box on
%     grid on
%     set(gca,'Layer','top')
%      filename = (['Psi3',num2str(i)]);
%          pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%     print(pt,'-dpng')
% 
%     figure(length(L0)+1)
%     axis off
%     colormap(crameri('nuuk',length(levels_cr)-1));
%     caxis([levels_cr(1),levels_cr(end)])
%     c=colorbar(gca,'west',TickLabelInterpreter='latex');
%     c.Label.Interpreter = 'latex';
%     c.Label.String      = '$log_{10}\left(\xi^{UM}\right), [n.d.]$';
% 
%     filename = (['Cbar1']);
%         pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%     print(filename,'-dpng')
% 
%     figure(length(L0)+2)
%     axis off
%     colormap(crameri('bilbao',length(level_la)-1));
%     caxis([level_la(1),level_la(end)])
%     c=colorbar(gca,'south',TickLabelInterpreter='latex');
%     c.Label.Interpreter = 'latex';
%     c.Label.String      = '$log_{10}\left(\eta_{eff}\right), [Pas]$';
%     filename = (['Cbar2']);
%         pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%     print(pt,'-dpng')
% 
%      figure(length(L0)+2)
%     axis off
%     colormap(crameri('berlin',length(level_Psi)-1));
%     caxis([level_Psi(1),level_Psi(end)])
%     c=colorbar(gca,'south',TickLabelInterpreter='latex');
%     c.Label.Interpreter = 'latex';
%     c.Label.String      = '$log_{10}\left(\frac{\Psi}{1+\xi^{UM}}\right), [Pas]$';
%     filename = (['Cbar3']);
%         pt=fullfile(folder_supplementary,filename);
%     set(gcf,'Color','w')
%     print(pt,'-dpng')
% 

%end

%%
