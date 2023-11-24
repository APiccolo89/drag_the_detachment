%=========================================================================%
% Realistic data: Script that run a few simulation based on the
% rheological data that have been compute in
% Compute_Viscosity_Olivine.m
%I have 4 parameters: Ts,Vd,Vn,L0; The other rheological parameter are
%given and considered constant together with the Pr,Tp.
%        eta0DS = f(Ts,Vd)
%        xiUS   = f(Ts,Vd,Vn,tau0)
%        eta0DM  = f(Tp,Vd,L0)
%        xiUM    = f(Tp,Vd,Vn,tau0)
%        tau0    = f(dRho,L0,g)
%        drho    = f(Ts,Tp)
%        L0      = K;
%        Vd/n      = K;
%        Ts        = K;
%        -> 1st The combination of these parameters gives situation
%        that are suitable and situation that are not solvable with the
%        current 0D approximation.
%        -> 2nd Certain parameter in conjuction with other gives
%        success or insuccess.
%        -> 3rd Q: What are the ranges of parameters that always gives
%        success?
%        -> 4th Q: What are the parameters that always gives a success,
%        gives reasonable td?
%        -> 5th Q: are these parameters reasonable?
%        .> 6th Q: From the point of view of a dislocation-diffusion
%        creep model, and on the insights provided by data, what are
%        the condition that allows a fast detachment (i.e. 10^-1-10^2)?
%=========================================================================%
clear all;
close all;
% add the path to the main function
addpath('Utilities/')
addpath('Utilities/Plot_Class/')
addpath('Utilities/ScientificColourMaps8/')
%load the Data Base
load('..\Data_Base\Realistic_Data_Base\Data_Base_Realistic_Scenario_1350_stress_L0_3D.mat');
% Main folder where to save the picture
folder = '../Manuscript_figure_folder/';
ptsave = '../Realistic_Cases/';
% General information of the size of each figure
if ~isfolder(ptsave)
    mkdir(ptsave);
end
% General information that are necessary (to change to a connected
% structure to the Data base
Tp    = S_initial_data.Tp;
Pr    = S_initial_data.Pr;
s0    = S_initial_data.s0;
if isnan(s0)
    flag = 1;
else
    flag=0;
end
Vnv   = S_initial_data.Vnv;
Vdv   = S_initial_data.Vdv;
T_mean   = S_initial_data.T_mean;
D0 = S_initial_data.D0;
L0 = S_initial_data.L0;
n=3.5;


%% ======================================================================= %
% Slab
%%=========================================================================%

% Create the vector to select part of the data
vnv = (2.0:1.0:27.0).*10^(-6);
vdv = (2.0:1.0:12.0).*10^(-6);
ind_n = zeros(length(vnv),1);
ind_d = zeros(length(vdv),1);
%find the index that correspond to the selected value
for i = 1:length(vnv)
    % Short comment: it appears that matlab has struggled to with
    % 15*10^(-6) and 15e-6. Their difference is above the numerical
    % limit, so, I needed to use this weird stuff to prevent any
    % problem. I simply give up on understanding.
    ind_n(i)=find(abs(Vnv-vnv(i))<1e-20,1);
end
for i = 1:length(vdv)
    ind_d(i)=find(Vdv==vdv(i));
end
% create a grid with all the data to facilitate the computation of the
% forward modelling
[ind_nm,ind_dm,L0m,T_meanm] = ndgrid(ind_n,ind_d,L0,T_mean);
success = ind_nm(:).*0.0;
% create additional matrix that contains some data
xium_      = ind_nm.*0.0;
eta0D_     = xium_;
% compute the total number of test that are going to computed
data_n = length(ind_nm(:));
disp(['Number of test is', num2str(data_n)])
n_success = 0;
for it = 1:data_n
    % Create the test name
    T_name   = strcat('T_',num2str(it));
    % select the index
    in = ind_nm(it);
    id = ind_dm(it);
    % Select the length
    l0 = L0m(it);
    % Select the temperature of the slab
    T_Slab = T_meanm(it);
    % Retrieve the viscosities of the slab and the non linear
    % coefficient
    if flag==1
        eta0DS = SLAB.eta0DS(in,id,L0 == l0,T_mean == T_Slab);
        xiUS   = SLAB.xiuS(in,id,L0 == l0,T_mean == T_Slab);
        eta0DM = UPPER_MANTLE.eta0DMP(in,id,L0 == l0,T_mean == T_Slab);
        xiUM   = UPPER_MANTLE.xiumP(in,id,L0 == l0,T_mean == T_Slab);
        s0     = UPPER_MANTLE.MMs0(in,id,L0 == l0,T_mean == T_Slab);
    else
        eta0DS = SLAB.eta0DS(in,id,T_mean == T_Slab);
        xiUS   = SLAB.xiuS(in,id,T_mean == T_Slab);
        eta0DM = UPPER_MANTLE.eta0DMP(in,id,L0 == l0);
        xiUM   = UPPER_MANTLE.xiumP(in,id,L0 == l0);
    end
    AVn(it) = Vnv(in);
    AVd(it) = Vdv(id);


    % Retrieve the average viscosity of the mantle as a function of the
    % lenght and the relative non linear coefficient
    tc(it)=nan;
    td(it)=nan;
    lambda(it) = nan;

    nlm = Problem_type;
    nlm.Linear=0;   % Switching the position of linear-non_linear activate the non linear upper mantle routine.
    nlm.iteration = 1;
    nlm.cut_off   = 0;
    % Call the function that process the simulation
    [Temp]   = Processing_simulation(eta0DS,xiUS,n,l0,s0,D0,eta0DM,0,xiUM,nlm);
    success(it)=0.0;
    % Show the progress
    % disp([num2str(it),'out of',num2str(data_n), 'tests'])
    % Select the the test that feature a initial Psi that is lower than
    % 1, and the not failed test
    if isfield(Temp,'D_norm')
        disp(Temp.initial_data.Psi)
        % Fill the test and the relative Metadata
        Tests.(T_name)=Temp;
        Meta_dataReal.Vd = Vdv(id);
        Meta_dataReal.Vn = Vnv(in);


        Meta_dataReal.T_Slab=T_Slab;
        Meta_dataReal.Pr    = Pr;
        Meta_dataReal.Tp    = Tp;
        [Meta_dataReal.Cd,Meta_dataReal.Cn,Meta_dataReal.phi] = UM.Compute_Cd_Cn(Pr,Vnv(in),Vdv(id),Tp);
        Meta_dataReal.w    = UM.rho.*UM.g;
        Tests.(T_name).Meta_dataReal = Meta_dataReal;
        % Disp the data of the success test
        disp(['Detachment_time = ',num2str(Tests.(T_name).t_det*3.5)]);
        disp(['T = ', num2str(T_Slab-273.15)])
        disp(['log10 etaM = ', num2str(log10(eta0DM))])
        disp(['log10 etaS = ', num2str(log10(eta0DS))])
        disp(['log10 xiUM = ', num2str(log10(xiUM))])
        disp(['log10 xiUS = ', num2str(log10(xiUS))])
        disp(['Vd         = ', num2str(Vdv(id)*1e6)])
        disp(['Vn = ', num2str((Vnv(in)*1e6))])
        disp(['l0 = ', num2str((l0./1e3))])
        success(it) = 1.0;
        tc(it)=Temp.initial_data.tc;
        td(it)=Temp.t_det;
        lambda(it) = Temp.Lambda(1);
    end
end
%%
% Structure containing all the combination of parameters:
% This structure contains the following parameters:
% a. The prime variable vectors that represents all the potential
% combination: i.e.
% vd x !Set of parameters that allow the generation of the fundamental
% vn y !variable of the 0D numerical experiments that is computable
% Ts z !   -> Vd-Vn grid sum up all the success. -> max success rate
% L0 u !   -> Vd-Vn Ts,L0 are free parameter that I,testing. max rate of
% S  1 !      success rate = L0*Ts
% =======================================================================%
Meta_data_Set_Tests = struct('Vnv',AVn,'Vdv',AVd,'TSlab',T_meanm,'L0',L0m, ...
    'max_success',length(L0)*length(T_mean),'Shape',[length(ind_n),length(ind_d),length(L0),length(T_mean)], ...
    'suc',success,'Lambda',lambda,'tc',tc,'td',td,'secMyear',365.25*60*60*24*1e6);
[Data_S] = extract_information_detachment(Tests,0,nlm);
save('../Data_Base/REAL_DATA_tp2','Tests','Data_S','Meta_data_Set_Tests')
%% 