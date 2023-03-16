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
addpath Realistic_Case/
addpath \Users\'Andrea Piccolo'\Dropbox\freezeColors-master\
close all
set(0,'defaultTextInterpreter','latex'); %trying to set the default
% Reference values:
Tp    = 1350+273.15;
Pr    = 3300*100e3*9.81;
t0    = 100e6;
Vnv   = [0e-6:0.1e-6:30e-6];
Vdv   = Vnv;
T_mean   = [600:25:1100]+273.15;
D0 = 80e3;
L0 = [300,400,500,600].*1e3;

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
[UPPER_MANTLE,SLAB] =  main_function_Real(t0, T_mean, Tp,Pr, D0,L0,Vnv,Vdv,UM,S);
% Wet Olivine Data:
[UPPER_MANTLE2,~] =  main_function_Real(t0, T_mean, Tp,Pr, D0,L0,Vnv,Vdv,UM2,S);

%% 
% Dry and Wet database 
save("Data_Base_Realistic_Scenario.mat", "UPPER_MANTLE2","SLAB","UPPER_MANTLE","UM2","UM");

%%





%%

for i = 1:length(L0)
    clf; close all;
    figure(i)

    xiuM = squeeze(UPPER_MANTLE2.xiumP(:,:,i));
    eta_eff_REF = (1./squeeze(UPPER_MANTLE2.eta0DMP(:,:,i))+1./squeeze(UPPER_MANTLE2.eta0MP(:,:,i))).^(-1);

    eta_eff_0  = squeeze(UPPER_MANTLE2.eta0DMP(:,:,i))./(1+xiuM.*eps(1).^((UM2.n-1)./UM2.n));
    eta_eff_5  = squeeze(UPPER_MANTLE2.eta0DMP(:,:,i))./(1+xiuM.*eps(2).^((UM2.n-1)./UM2.n));

    levels_cr = 1:1:10;
    level_la  = 17:1:26;

    tiledlayout(1,4, 'TileSpacing','compact', 'Padding', 'none');
    nexttile;
    a=pcolor(x,y,log10(xiuM));shading interp;
    colormap(crameri('nuuk',length(levels_cr)-1));
    xlabel('$V_d [10^6\frac{m^3}{J}]$',Interpreter='latex')
    ylabel('$V_n [10^6\frac{m^3}{J}]$',Interpreter='latex')
    axis_x= get(gca, 'XAxis');
    axis_x.TickLabelInterpreter = 'latex';
    axis_x.TickValues   = [5:5:25];

    axis_y= get(gca, 'YAxis');
    axis_y.TickLabelInterpreter = 'latex';
    axis_y.TickValues   = [5:5:25];


    title(['$\xi^{UM}, L_0 = ',num2str(L0(i)./1000),' [km]$'],Interpreter='latex')
    caxis([levels_cr(1),levels_cr(end)])
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
    line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
    line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
    hold off
    axis square;
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    nexttile

    c=pcolor(x,y,log10(eta_eff_0));shading interp;
    colormap(crameri('bilbao',length(level_la)-1));
    xlabel('$V_d [10^6\frac{m^3}{J}]$',Interpreter='latex')
    title('$v = 0 [\frac{cm}{yrs}]$',Interpreter='latex')
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
    axis square;
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;

    nexttile
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
    title('$v = 5 [\frac{cm}{yrs}]$',Interpreter='latex')
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
    line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
    line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
    hold off
    caxis([level_la(1),level_la(end)])
    axis square;
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;

    nexttile
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
    title('Reference condition',Interpreter='latex')
    hold on
    line([2,2],[2,27],'Color' ,'r',LineWidth = 1.0)
    line([10,10],[2,27],'Color', 'r',LineWidth = 1.0)
    line([2,10],[2,2],'Color' ,'r',LineWidth = 1.0)
    line([2,10],[27,27], 'Color' ,'r',LineWidth = 1.0)
    hold off
    caxis([level_la(1),level_la(end)])
    axis square;
    box on
    grid on
    set(gca,'Layer','top')
    freezeColors;
    filename = (['Wet_Mantle',num2str(L0(i)./1e3),'.png']);
    print(filename,'-dpng')

    figure(length(L0)+1)
    axis off
    colormap(crameri('nuuk',length(levels_cr)-1));
    caxis([levels_cr(1),levels_cr(end)])
    c=colorbar(gca,"west",TickLabelInterpreter="latex");
    c.Label.Interpreter = 'latex';
    c.Label.String      = '$log_{10}\left(\xi^{UM}\right), [n.d.]$';

    filename = (['Cbar1']);
    print(filename,'-dpng')

    figure(length(L0)+2)
    axis off
    colormap(crameri('bilbao',length(level_la)-1));
    caxis([level_la(1),level_la(end)])
    c=colorbar(gca,"south",TickLabelInterpreter="latex");
    c.Label.Interpreter = 'latex';
    c.Label.String      = '$log_{10}\left(\eta_{eff}\right), [Pas]$';
    filename = (['Cbar2']);
    print(filename,'-dpng')


end

%%
clf; close all;
figure(5)

    axis off

   z_min = min(Vnv)*1e6;  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = max(Vnv)*1e6;  %log10(max(Data_S(1,:)));
    % See if the user installed Crameri cmap utilities, otherwise punish
    % him with jet colormap by default
    % Shamelessly copied from 
    % https://de.mathworks.com/matlabcentral/answers/595864-how-to-plot-a-line-graph-x-y-with-a-color-bar-representing-z-axis
    
    try
        cmap = colormap(crameri('nuuk'));
    catch
        cmap = colormap('jet');
    end
    % Set colorbar
    c=colorbar(gca,"east",TickLabelInterpreter="latex");
    c.Label.String = '$V_{n|d} [10^6 \frac{m^3}{J}]$';
    c.Label.Interpreter = 'latex';
    % Define coloraxis
    caxis([z_min z_max])
    % min/max are rounded using the log values of Lambda
    val = (z_max-z_min)/10;
    val2 = z_min:val:z_max;

    Ticks=val2;
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('{%0.4g}',x),Ticks,'UniformOutput',false);
    try
        set(c,'Ticks',Ticks);
        set(c,'TickLabels',TickLabels);
    catch %HG1
        TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
        set(c,'YTick',Ticks);
        set(c,'YTickLabel',TickLabels);
    end
    filename = (['SlabCbar']);
print(filename,'-dpng')

    % Set colorbar
   figure(6)
    for i = 1:length(Vnv)
 
         V = (Vnv(i)*1e6-z_min)/(z_max-z_min);
        if V<0 
              V=0;
        elseif V>1
          V=1;
        end
        V=round(1+V*(size(cmap,1)-1));%round to nearest index
        C = cmap(V,:);
        hold on
        plot(squeeze(SLAB.MMT_av(1,1,:))-273.15,(squeeze(SLAB.eta0S(1,i,:))),'Color',C,'LineStyle','-','LineWidth',0.8);
        plot(squeeze(SLAB.MMT_av(1,1,:))-273.15,(squeeze(SLAB.eta0DS(i,1,:))),'Color',C,'LineStyle','-','LineWidth',0.8);
        set(gca, 'YScale', 'log')
        xlabel('$T [^\circ C]$',Interpreter='latex')
        ylabel('$log_{10}\left(\eta^{S}\right) [Pas]$',Interpreter='latex')

    end
        plot(squeeze(SLAB.MMT_av(1,1,:))-273.15,(squeeze(SLAB.eta0S(1,end,:))),'Color','k','LineStyle','--','LineWidth',2.0);
        plot(squeeze(SLAB.MMT_av(1,1,:))-273.15,(squeeze(SLAB.eta0S(1,1,:))),'Color','k','LineStyle','--','LineWidth',2.0);

box on
grid on
plot(squeeze(SLAB.MMT_av(1,1,:))-273.15,(squeeze(SLAB.eta0DS(end,10:10:end,:))),'Color','r','LineStyle','--','LineWidth',2.0);
plot(squeeze(SLAB.MMT_av(1,1,:))-273.15,(squeeze(SLAB.eta0DS(1,10:10:end,:))),'Color','r','LineStyle','--','LineWidth',2.0);
a = (squeeze(SLAB.MMT_av(1,1,:))./squeeze(SLAB.MMT_av(1,1,:)))*1e24;
b = a./1e6;
x= squeeze(SLAB.MMT_av(1,1,:));
plot(x-273.15,a,'Color','k',LineWidth=2.0)
plot(x-273.15,b,'Color','k',LineWidth=2.0)
fill([x-273.15,x-273.15],[a,fliplr(b)],'r',FaceAlpha=0.5)
set(gca, 'YScale', 'log')
axis square
filename = (['Slab']);
print(filename,'-dpng')


clf; close all;
%%



eta_n = (1./SLAB.eta0S+1./SLAB.eta0DS).^(-1);
t_c   = 2.*eta_n./t0;
t_c   = t_c./(365.25.*24.*60.*60.*3.5);
 ind = find(Vnv >=2.0e-6 & Vnv <=27e-6);
ind2 = find(Vdv >=2.0e-6 & Vdv <=12e-6);
t_cA = t_c(ind,ind2,:);
t_cB = zeros(length(T_mean),1);
t_cM = t_cB; 
t_cm = t_cB; 
eta_M = t_cB; 
eta_m = t_cB;
eta_av = t_cB; 
for i = 1:length(T_mean)
    B=squeeze(eta_n(:,:,i));
    eta_av(i) = mean(B(:)); 
    eta_m(i) = min(B(:)); 
    eta_M(i) =max(B(:)); 
    t_cB(i)= (2.*eta_av(i)./t0)./(365.25*60*60*24);
        t_cM(i)= (2.*eta_M(i)./t0)./(365.25*60*60*24);
    t_cm(i)= (2.*eta_m(i)./t0)./(365.25*60*60*24);

end

figure(7)
yyaxis left

plot(T_mean-273.15,t_cB./1e6,LineStyle='-',LineWidth=2.0);
hold on
plot(T_mean-273.15,t_cm./1e6,LineStyle='-',LineWidth=0.8);
plot(T_mean-273.15,t_cM./1e6,LineStyle='-',LineWidth=0.8);
box on
grid on
set(gca, 'YScale', 'log')
%ylim([0,100])
ylabel('$t_d [Myrs]$', Interpreter='latex')
xlabel('$T [^\circ C]$', Interpreter='latex')
yyaxis right
plot(T_mean-273.15,eta_av,LineStyle='-',LineWidth=2.0);
hold on
plot(T_mean-273.15,eta_m,LineStyle='-',LineWidth=0.8);
plot(T_mean-273.15,eta_M,LineStyle='-',LineWidth=0.8);
axis square
set(gca, 'YScale', 'log')

filename = (['Slab_tc']);
print(filename,'-dpng')

%%
ind = find(Vnv >=2.0e-6 & Vnv <=27e-6);
ind2 = find(Vdv >=2.0e-6 & Vdv <=12e-6);
xiuM_choosen = (UPPER_MANTLE.xiumP(ind2,ind,:));
eta0D_choosen = (UPPER_MANTLE.eta0DMP(ind2,ind,:));

xiuS_choosen = (SLAB.xiuS(ind2,ind,:));
eta0DS_choosen = (SLAB.eta0DS(ind2,ind,:));
eta0S_choosen = (SLAB.eta0S(ind2,ind,:));
t_c_choosen   = (t_c(ind2,ind,:))

Lambda  = zeros(length(ind2),length(ind),length(L0),length(T_mean));
Psi_0   = Lambda; 
for i =1:length(L0)
    for j = 1:length(T_mean)
    l = L0(i); 
    gamma  = (l*5)/(2000e3);
    etan = ((1./eta0S_choosen(:,:,j)+1./eta0DS_choosen(:,:,j))).^-1;
    Psi_0(:,:,i,j)    = squeeze(eta0D_choosen(:,:,i))./squeeze(etan(:,:,1)); 
    Lambda(:,:,i,j) = (Psi_0(:,:,i,j).*gamma)./(1+squeeze(xiuM_choosen(:,:,i)));
    end
end


vla = 0.0;
figure(11)
clf;
pd = fitdist(log10(Lambda(:)),'normal');
% Find the pdf that spans the disribution
x_pdf = linspace(min(log10(Lambda(:))),max(log10(Lambda(:))));
y_pdf = pdf(pd,x_pdf);

h=histogram(log10(Lambda(:)),50,'Normalization','pdf');
A = squeeze(Lambda(:,:,1,:));
B = squeeze(Lambda(:,:,4,:));
hold on

h.FaceAlpha = 0.2; 
xlabel('$log_{10}\left(\frac{\Lambda}{1+\xi^{UM}}\right), [n.d.]$',Interpreter='latex')
ylabel('$f, [n.d.]$',Interpreter='latex')

grid on;
box on;
filename = (['Hist']);
print(filename,'-dpng')



%% Figure Upper Mantle 
% Prepare data: 
% Lambda0 = (1+xiS)eta0DM/(1+xium)eta0DS*gamma 
%
gamma  = (L0*5)/(2000);
T_plot   = 900; 
t_c      = SLAB.eta0DS(10,:)./100e6;

F1=figure(1)
% plot xium vs L0 @ constant vd 
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
 
    

% compute Psi_0 from a slab that has fixed average temperature 
iT = find(T_mean==T_plot); 
xiuS = squeeze(SLAB.xiuS(:,:,iT));
eta0DS = squeeze(SLAB.eta0DS(:,:,iT));
eta0S  = squeeze(SLAB.eta0S(:,:,iT));
etan = 1./(1./eta0DS+1./eta0S);
eta0DM = squeeze(UPPER_MANTLE.eta0DMP(1,:,:))./(1+squeeze(UPPER_MANTLE.xiumP(1,:,:)));
Lambda_0 = eta0DM.*0;
Psi_0    = Lambda_0;
for i = 1:1:length(L0)
    Psi_0(:,i)    = eta0DM(:,i)./etan'; 
    Lambda_0(:,i) = Psi_0(:,i).*gamma(i); 
end


F2=figure(2)
% plot xium vs L0 @ constant vd 
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = min(Vnv);  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = max(Vnv);  %log10(max(Data_S(1,:)));
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
   
    for i = 1:length(Vnv)
         c=colorbar;
    c.Label.String = 'V_a []';
    % Define coloraxis
    caxis([z_min z_max])
    % min/max are rounded using the log values of Lambda
    val = (z_max-z_min)/10;
    val2 = z_min:val:z_max;

    Ticks=val2;
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('{%0.4g}',x),Ticks,'UniformOutput',false);
    try
        set(c,'Ticks',Ticks);
        set(c,'TickLabels',TickLabels);
    catch %HG1
        TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
        set(c,'YTick',Ticks);
        set(c,'YTickLabel',TickLabels);
    end
         V = (Vnv(i)-z_min)/(z_max-z_min);
        if V<0 
              V=0;
        elseif V>1
          V=1;
        end
        V=round(1+V*(size(cmap,1)-1));%round to nearest index
        C = cmap(V,:);
        hold on
        plot(L0/1000, Lambda_0(i,:),"Color",C)
        set(gca, 'YScale', 'log')
        xlabel('$L_0 [km]$',Interpreter='latex')
        ylabel('$\Lambda_{0} [n.d.]$',Interpreter='latex')
        title('$\Lambda_{0} [n.d.], T_P = 1250 [^\circ C], Vd = 5e-6,T=1200 K$',Interpreter='latex')

    end
    hold off
    box on
    grid on
F3=figure(3)
% plot xium vs L0 @ constant vd 
    % Collect the field names 
    % Set the min and max of lambda value for the coloring of the plot
    z_min = min(Vnv);  %round(log10(min(Data_S(1,Data_S(1,:)<1.0))));%log10(min(Data_S(1,:)));
    z_max = max(Vnv);  %log10(max(Data_S(1,:)));
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
   
    for i = 1:length(Vnv)
         c=colorbar;
    c.Label.String = 'V_a []';
    % Define coloraxis
    caxis([z_min z_max])
    % min/max are rounded using the log values of Lambda
    val = (z_max-z_min)/10;
    val2 = z_min:val:z_max;

    Ticks=val2;
    % Produce the tick
    TickLabels=arrayfun(@(x) sprintf('{%0.4g}',x),Ticks,'UniformOutput',false);
    try
        set(c,'Ticks',Ticks);
        set(c,'TickLabels',TickLabels);
    catch %HG1
        TickLabels=strrep(strrep(TickLabels,'{',''),'}','');%remove TeX formatting
        set(c,'YTick',Ticks);
        set(c,'YTickLabel',TickLabels);
    end
         V = (Vnv(i)-z_min)/(z_max-z_min);
        if V<0 
              V=0;
        elseif V>1
          V=1;
        end
        V=round(1+V*(size(cmap,1)-1));%round to nearest index
        C = cmap(V,:);
        hold on
        plot(L0/1000, Psi_0(i,:),"Color",C)
        set(gca, 'YScale', 'log')
        xlabel('$L_0 [km]$',Interpreter='latex')
        ylabel('$\Psi_{0} [n.d.]$',Interpreter='latex')
        title('$\Psi_{0} [n.d.], T_P = 1250 [^\circ C], Vd = 5e-6,T=1200 K$',Interpreter='latex')

    end
    hold off
    box on
    grid on


