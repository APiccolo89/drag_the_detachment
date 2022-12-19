function [P_Var,D,t,Tau,d,tdet] = Reading_Data_Base(Test_name,DB_path,Df_S)
%==========================================================================
% Function to read h5 file data.
% Input : Test Name, Data base path
% Output: Primary Variables {i.e. the variable that can be accepted by
%         compute the slab characteristics}
%==========================================================================
% Retrieve the information for replicate the 2D tests and to replicate them
% using the dimensionless 0D equation
%==========================================================================
% Path Construction
status = strcat(Test_name,'/failed');

P_Var.failed = h5read(DB_path,status);

path_D = strcat(Test_name,'/TimeEvolution/Slab1D/D');
path_t = strcat(Test_name,'/TimeEvolution/time');
path_L0 = strcat(Test_name,'/Initial_Data/Initial_Condition/L0');
path_tc = strcat(Test_name,'/Initial_Data/Initial_Condition/tc');
path_s0 = strcat(Test_name,'/Initial_Data/Initial_Condition/tau0');
path_eta0DS = strcat(Test_name,'/Initial_Data/Initial_Condition/eta0');
path_eta0DM = strcat(Test_name,'/Initial_Data/Initial_Condition/Astenosphere/eta');
path_tzz = strcat(Test_name,'/TimeEvolution/Slab1D/tauzz');
path_txx = strcat(Test_name,'/TimeEvolution/Slab1D/tauxx');
path_txz = strcat(Test_name,'/TimeEvolution/Slab1D/tauxz');%/Viscous/High_Resolution_NLM/T_I0_VA20_V22_LM1_NE/Detachment/det_vec/
path_depth = strcat(Test_name,'/Detachment/det_vec/Depth');
path_time_det = strcat(Test_name,'/Detachment/det_vec/t_det_td');
P_Var.L0 = h5read(DB_path,path_L0)*1e3;
P_Var.tc =  h5read(DB_path,path_tc);
P_Var.s0 =  h5read(DB_path,path_s0)*1e6;
P_Var.eta0DS = Df_S*h5read(DB_path,path_eta0DS);
P_Var.eta0DM = h5read(DB_path,path_eta0DM);


disp('=====================================================================')
disp(['node ',Test_name, 'is processed'])
if P_Var.failed ==0
    disp(['it was a successful test'])
else
    disp(['It was not a succesful test'])
end
disp('=====================================================================')
disp(['L0 is ', num2str(P_Var.L0./1e3), ' [km]'])
disp(['tc is ', num2str(P_Var.tc), ' [Myrs]'])
disp(['tau0 is ', num2str(P_Var.s0./1e6), ' [MPa]'])
disp(['log10(eta0DS) is ', num2str(log10(P_Var.eta0DS)), ' [Pa.s]'])
disp(['log10(eta0DM) is ', num2str(log10(P_Var.eta0DM)), ' [Pa.s]'])
disp('=====================================================================')
disp('=====================================================================')
disp('When a test is failed, means that the detachment scales are higher than 10 n*t/tc')
disp('or that the detachment timescale are exceeding 100 Myrs, which has no geological sense')


if P_Var.failed == 1
    D = [];
    t= [];
    Tau = [];
    d  = [];
    tdet = []; 
else
    tdet   =  h5read(DB_path,path_time_det);

    % Relevant Information
    D = h5read(DB_path,path_D);
    D=D./D(1);
    t= h5read(DB_path,path_t);
    t = t./P_Var.tc;
    tau_zz = h5read(DB_path,path_tzz);
    tau_xx = h5read(DB_path,path_txx); 
    tau_xz = h5read(DB_path,path_txz);
    Tau = [tau_zz(2),tau_xx(2),tau_xz(2)]./(P_Var.s0./1e6);
    d = (h5read(DB_path,path_depth)+100)./(P_Var.L0/1000); 
        
end

end
