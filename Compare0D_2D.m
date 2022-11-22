%Function 2 read the data base 
%Extract data 
%Extract primary variable (i.e. L0, eta0DM,eta0DS,tau0)
%call slab characteristic -> run 0D (fetch parmater x->x*L0 )

clear all; 
close all;
addpath Adimensional\
addpath Dimensional\
addpath Utilities\

clear all 
% Common Data
D0=80e3; 
n=3.5;
Df_S=10; 
nlm = Problem_type.Linear;
Df_UM = nan; 

DB_path = '../Test_DB.hdf5';
l = h5info(DB_path,'/Viscous/High_Resolution_NLM');
Tests  = {l.Groups.Name};
Chosen = Tests{10};
[P_Var,D,t] = Reading_Data_Base(Chosen,DB_path,Df_S);




[ID] = Compute_slab_characteristics(P_Var.eta0DS,Df_S,n,P_Var.L0,P_Var.s0,D0,P_Var.eta0DM,Df_UM,nlm);
% Example
f = 210.011;
ID.ID_A.Lambda = ID.ID_A.Lambda*f;
[TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
time_0D = TestData.time; 
D_0D    = TestData.D_norm;
D_0D2D = interp1(time_0D,D_0D,t,'linear','extrap');
D_0D2D(isnan(D_0D2D) | D_0D2D<0.1)=0.0;
D(isnan(D) | D<0.1)=0.0;
fit = goodnessOfFit(D,D_0D2D,'NRMSE');
disp(['fit ',num2str(fit*100),'%'])

figure(1)
plot(n.*t,D,'Color','k',LineStyle='-.',LineWidth=1.2)
hold on
plot(n.*t,D_0D2D,'Color','b',LineStyle='-.',LineWidth=1.2)
xlim([0,10])
ylim([0.1,1.0])
legend('2D numerical experiment','0D numerical experiment')
grid on 
xlabel('$\frac{n\cdot t}{t_c}$',Interpreter='latex')
ylabel('$\frac{D}{D_0}$',Interpreter='latex')
figure(2)
plot(n.*t,D-D_0D2D)
xlim([0,10])
grid on 
xlabel('$\frac{n\cdot t}{t_c}$',Interpreter='latex')
ylabel('$\frac{D_{2D}-D_{0D}}{D_0}$',Interpreter='latex')
function [P_Var,D,t] = Reading_Data_Base(Test_name,DB_path,Df_S)
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
path_D = strcat(Test_name,'/TimeEvolution/Slab1D/D');
path_t = strcat(Test_name,'/TimeEvolution/time');
path_L0 = strcat(Test_name,'/Initial_Data/Initial_Condition/L0');
path_tc = strcat(Test_name,'/Initial_Data/Initial_Condition/tc');
path_s0 = strcat(Test_name,'/Initial_Data/Initial_Condition/tau0');
path_eta0DS = strcat(Test_name,'/Initial_Data/Initial_Condition/eta0');
path_eta0DM = strcat(Test_name,'/Initial_Data/Initial_Condition/Astenosphere/eta');
% Relevant Information
D = h5read(DB_path,path_D);
D=D./D(1);
t= h5read(DB_path,path_t); 
P_Var.L0 = h5read(DB_path,path_L0)*1e3; 
P_Var.tc =  h5read(DB_path,path_tc);
P_Var.s0 =  h5read(DB_path,path_s0)*1e6;
P_Var.eta0DS = Df_S*h5read(DB_path,path_eta0DS);
P_Var.eta0DM = h5read(DB_path,path_eta0DM);
P_Var.failed = h5read(DB_path,status);
t = t./P_Var.tc; 
disp('=====================================================================')
disp(['node ',Test_name, 'is processed'])
if P_Var.failed ==0
disp(['it was a successful test'])
else
disp(['It was not a succesful test'])
end
disp('=====================================================================')
end
