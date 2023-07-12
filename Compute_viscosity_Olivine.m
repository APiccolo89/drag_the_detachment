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

% Dry and Wet database 
save((filename), 'UPPER_MANTLE2','SLAB','UPPER_MANTLE','UM2','UM','S_initial_data');



