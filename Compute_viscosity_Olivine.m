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

% Reference values: 
Tp    = 1250+273.15; 
Pr    = 3300*100e3*9.81; 
t0    = 50e6;
Vnv   = [0e-6:1e-6:20e-6];
Vdv   = 5e-6; 
T_mean   = [700:100:1300]; 
D0 = 80e3;
L0 = [1:10:600].*1e3;
%Cn = (En*phi-Vnv*(1-phi*Pr))/(R*Tp);
%Cd = (Ed*phi-Vd*(1-phi*Pr))/(R*Tp);
UM = Mantle_Unit_Properties(3300,3e-5,1050,1.5e-9,6.22254e-16,375e3,530e3,3.5);
S         = Mantle_Unit_Properties(3360,3e-5,1050,1.5e-9,6.22254e-16,375e3,530e3,3.5);
% Dry Olivine Data: 
[UPPER_MANTLE,SLAB] =  main_function_Real(t0, T_mean, Tp,Pr, D0,L0,Vnv,Vdv,UM,S);
l10xium = squeeze(UPPER_MANTLE.xiumP(1,:,:));
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
        plot(L0/1000, l10xium(i,:),"Color",C)
        set(gca, 'YScale', 'log')
        xlabel('$L_0 [km]$',Interpreter='latex')
        ylabel('$\xi^{UM} [n.d.]$',Interpreter='latex')
        title('$\bar{\xi}^{UM}, T_P = 1250 [^\circ C], Vd = 5e-6$',Interpreter='latex')

    end
    hold off
    box on
    grid on
    

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


