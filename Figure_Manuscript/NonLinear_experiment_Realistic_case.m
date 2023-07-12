%=========================================================================%
%  Manuscript Non Linear Experiments 
%=========================================================================%

close all;
clear all;
clf;
addpath ../Realistic_Case
addpath ../Utilities
addpath \Users\'Andrea Piccolo'\Dropbox\freezeColors-master\
addpath \Users\Andrea Piccolo\Dropbox\Bayreuth_Marcel\export_fig
addpath C:\Users\AndreaPiccolo\Downloads\github_repo(1)
tic;
load('../../Data_Base/REAL_DATA_tp2');
toc;
pt_save=('../../Manuscript/Realistic_Cases');if ~isfolder(pt_save); mkdir(pt_save); end
%% Reshape the datainto matrices
Vn_test = reshape(Meta_data_Set_Tests.Vnv,Meta_data_Set_Tests.Shape);
Vd_test = reshape(Meta_data_Set_Tests.Vdv,Meta_data_Set_Tests.Shape);
TS_test = reshape(Meta_data_Set_Tests.TSlab,Meta_data_Set_Tests.Shape);
L0_tests      = reshape(Meta_data_Set_Tests.L0,Meta_data_Set_Tests.Shape);
suc_test     = reshape(Meta_data_Set_Tests.suc,Meta_data_Set_Tests.Shape);
td_test     = reshape(Meta_data_Set_Tests.td,Meta_data_Set_Tests.Shape);
Lambda_test     = reshape(Meta_data_Set_Tests.Lambda,Meta_data_Set_Tests.Shape);
tc_test     = reshape(Meta_data_Set_Tests.tc,Meta_data_Set_Tests.Shape);


font_axes = 16; 
font_legend = 14; 
font_text   = 5; 
size_picture = [12,12.5];
LineWidth = 1.0; 
marker_size = 10;

% Figure
figure(1)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])
clf; 
z =(Data_S.Lambda)./(1+Data_S.xiUM);  
size = 5+(40-5)./((max(z)-min(z))).*(z-min(z));
a = Data_S.tc./365.25./24./60./60./1e6./Data_S.n; 
b = ((Data_S.tdet)).*Data_S.tc./Data_S.n;
b = b./365.25./24./60./60./1e6; 
ax = gca; 
scatter(ax,a,Data_S.tdet,10,log10(z),'filled','MarkerEdgeColor','k');
colormap(crameri('bilbao',20));
ax.XScale = 'log';
ax.YScale = 'log'; 
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.LineWidth = 1.2;
ax.XLim = [10^-1,10^2];
caxis([-5 -1])
%ax.XTickLabel = []; 
%ax.YTickLabel = [];
%yline(100,'LineStyle','-','LineWidth',1.5);
namefigure='Figure_timescales.png';
pt = fullfile(pt_save,namefigure);
print(pt,'-dpng');



%%
allowed=b>10^(-1) & b<10^2; %& log10(Data_S.xiUM)>1.0;

figure(2)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])

clf; 
z =Data_S.tdet;  
size = 5+(40-5)./((max(z)-min(z))).*(z-min(z));
y = (Data_S.tau_max.*Data_S.tau0)./1e9; 
x = Data_S.Lambda./(1+Data_S.xiUM);
z = b;
l = Data_S.L0;
x = x(allowed);
y = y(allowed);
z = z(allowed);
l = l(allowed)./1000;
size = size(allowed);
ax = gca; 
hold on
p1=scatter(x(l==300),y(l==300),10,[176 	23 	31]./255,'filled','o','MarkerEdgeColor','k');
p2=scatter(x(l==400),y(l==400),10,[139 	71 	137]./255','filled','o','MarkerEdgeColor','k');
p3=scatter(x(l==500),y(l==500),10,[0 191 255]./255,'filled','o','MarkerEdgeColor','k');
p4=scatter(x(l==600),y(l==600),10,[0 0 0 ],'filled','o','MarkerEdgeColor','k');
l=legend('$L_0$=300, [km]','$L_0$=400, [km]','$L_0$=500, [km]','$L_0$=600, [km]');
l.Interpreter = 'Latex';
l.FontSize = 12;
l.Location = 'southwest';

colormap(crameri('oslo',10));
ax.XScale = 'log';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XLim = ([10^-5,10^0]);
caxis([-1 2])
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.LineWidth = 1.2;
ax.XTickLabel = []; 
ax.YTickLabel = [];
hold on 
%yline(1.8,'LineWidth',1.2,'Alpha',0.3)
bla =0; 
namefigure='Stress.png';
pt = fullfile(pt_save,namefigure);
print(pt,'-dpng');
% figure(2)
% colorbar(ax);
% delete(p1)
% namefigure='StressBar.png';
% ax.Visible = 'off';
% pt = fullfile(pt_save,namefigure);
% print(pt,'-dpng');
%%
figure(3)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])

clf; 
z =Data_S.tdet;  
size = 5+(40-5)./((max(z)-min(z))).*(z-min(z));
y = (Data_S.tau_max); 
x = Data_S.Lambda./(1+Data_S.xiUM);
z = b;
l = Data_S.L0;
x = x(allowed);
y = y(allowed);
z = z(allowed);
l = l(allowed)./1000;
size = size(allowed);
ax = gca; 
hold on
p1=scatter(x(l==300),y(l==300),10,[176 	23 	31]./255,'filled','o','MarkerEdgeColor','k');
p2=scatter(x(l==400),y(l==400),10,[139 	71 	137]./255','filled','o','MarkerEdgeColor','k');
p3=scatter(x(l==500),y(l==500),10,[0 191 255]./255,'filled','o','MarkerEdgeColor','k');
p4=scatter(x(l==600),y(l==600),10,[0 0 0 ],'filled','o','MarkerEdgeColor','k');
colormap(crameri('oslo',10));
ax.XScale = 'log';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XLim = ([10^-5,10^0]);
caxis([-1 2])
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.LineWidth = 1.2;
ax.XTickLabel = []; 
ax.YTickLabel = [];
hold on 
%yline(1.8,'LineWidth',1.2,'Alpha',0.3)
bla =0; 
namefigure='Stress_ND.png';
pt = fullfile(pt_save,namefigure);
print(pt,'-dpng');

%% Create the folder to save
% Prepare the 3D visualisation 
time_d_dim  = td_test.*tc_test./Meta_data_Set_Tests.secMyear; 
allowed     = suc_test(:) == 1 & time_d_dim(:)>0.1 & time_d_dim(:) <120.0; 
not_        = allowed==0 & ~isnan(Lambda_test(:)); 


figure(3)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])

clf; 
ax = gca; 
p1=scatter3(Vn_test(allowed),Vd_test(allowed),TS_test(allowed)-273.15,30,log10(Lambda_test(allowed)),"filled",'MarkerEdgeColor','k',"MarkerFaceAlpha",0.7); 
hold on 
p2 = scatter3(Vn_test(not_),Vd_test(not_),TS_test(not_)-273.15,6,log10(Lambda_test(not_)),"filled",'d',"MarkerFaceAlpha",0.2,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1); 
ax.XGrid = 'on'; 
ax.YGrid = 'on';
ax.ZGrid = 'on'; 
ax.LineWidth = 1.2; 
ax.Box = 'on';
ax.XTickLabel = []; 
ax.YTickLabel = [];
ax.ZTickLabel = [];
colormap(ax,crameri('bilbao',10)); 

%ax.BoxStyle = 'full';
view([-45 60 30])
caxis([-5,-1]);
fname = ['3D_TS_Plot.png'];
pt = fullfile(pt_save,fname);
print(pt,'-dpng')

figure(3)
ax.Visible = 'off';
colorbar(ax);
delete(p1);delete(p2);
fname = ['CBar_1.png'];
pt = fullfile(pt_save,fname);
print(pt,'-dpng')


figure(4)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])

clf; 
ax = gca; 
p1=scatter3(Vn_test(allowed),Vd_test(allowed),L0_tests(allowed)./1e3,30,log10(Lambda_test(allowed)),"filled",'MarkerEdgeColor','k',"MarkerFaceAlpha",0.7); 
hold on 
p2 = scatter3(Vn_test(not_),Vd_test(not_),L0_tests(not_)./1e3,6,log10(Lambda_test(not_)),"filled",'d',"MarkerFaceAlpha",0.2,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1); 
colormap(crameri('bilbao',10)); 
ax.XGrid = 'on'; 
ax.YGrid = 'on';
ax.ZGrid = 'on'; 
ax.LineWidth = 1.2; 
ax.XTickLabel = []; 
ax.YTickLabel = [];
ax.ZTickLabel = [];
ax.Box = 'on';
%ax.BoxStyle = 'full';
view([-45 60 30])
caxis([-5,-1]);
fname = ['3D_L0_Plot.png'];
pt = fullfile(pt_save,fname);
print(pt,'-dpng')

figure(5)
set(gcf, 'Units','centimeters', 'Position', [0, 0, size_picture(1),size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [size_picture(1), size_picture(2)])

clf; 
ax = gca; 
p1=scatter(Vn_test(allowed),TS_test(allowed)-273.15,30,log10(Lambda_test(allowed)),"filled",'MarkerEdgeColor','k'); 
hold on 
p2 = scatter(Vn_test(not_),TS_test(not_)-273.15,5,log10(Lambda_test(not_)),"filled",'d',"MarkerFaceAlpha",0.1,'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1); 
colormap(crameri('bilbao',10)); 
ax.Box = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on'; 
ax.XTickLabel = []; 
ax.YTickLabel = [];
    xline(15e-6)
xline(20e-6)
ax.LineWidth = 1.2; 
caxis([-5,-1]);
fname = ['2D_TS_Plot.png'];
pt = fullfile(pt_save,fname);
print(pt,'-dpng')

allowed2 = suc_test(:) == 1 & time_d_dim(:)>0.1 & time_d_dim(:) <120.0 & Vn_test(:) >=15e-6 & Vn_test(:) <=20e-6; 




