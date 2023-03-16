function [P_Var,D,t,Tau,d,tdet,topo,tauL] = Reading_Data_Base(Test_name,DB_path,Df_S,nlm)
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
if nlm.islinear == 0
    path_eta0M = strcat(Test_name,'/Initial_Data/Initial_Condition/Astenosphere/eta0');
end
path_eta0DM = strcat(Test_name,'/Initial_Data/Initial_Condition/Astenosphere/eta');
path_tau = strcat(Test_name,'/TimeEvolution/Slab1D/tau');
path_depth = strcat(Test_name,'/Detachment/det_vec/Depth');
path_time_det = strcat(Test_name,'/Detachment/det_vec/t_det_td');
path_topo     = strcat(Test_name,'/TimeEvolution/Free_Surface_Time/Amplitude');
path_tauL     = strcat(Test_name,'/TimeEvolution/Whole_Lithosphere/tauii');

P_Var.L0 = h5read(DB_path,path_L0)*1e3;
P_Var.tc =  h5read(DB_path,path_tc);
P_Var.s0 =  h5read(DB_path,path_s0)*1e6;
P_Var.eta0DS = Df_S*h5read(DB_path,path_eta0DS);
P_Var.eta0DM = h5read(DB_path,path_eta0DM);
P_Var.xiUM   =NaN;
if nlm.islinear == 0
    P_Var.eta0M = h5read(DB_path,path_eta0M);
    P_Var.xiUM   = P_Var.eta0DM./P_Var.eta0M; 
    disp(['xium is ', num2str(P_Var.xiUM)])
end


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
    topo = []; 
    tauL= []; 
else
    tdet   =  h5read(DB_path,path_time_det);

    % Relevant Information
    D = h5read(DB_path,path_D);
    D=D./D(1);
    t= h5read(DB_path,path_t);
    t = t./P_Var.tc;
    tau = h5read(DB_path,path_tau);

    Tau = tau./(P_Var.s0./1e6);
    d = (h5read(DB_path,path_depth)+100)./(P_Var.L0/1000); 
    topo_1  = h5read(DB_path,path_topo);
    topo    = mean(topo_1,2);
    tauL    = h5read(DB_path,path_tauL)./(P_Var.s0./1e6);
        
end

end
