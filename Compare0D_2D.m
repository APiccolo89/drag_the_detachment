%=========================================================================%

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
pt_save = '../Results/Comparison2D_0D';
if not(isfolder(pt_save))
     mkdir(pt_save);
end
l = h5info(DB_path,'/Viscous/High_Resolution_NLM');
Tests  = {l.Groups.Name};
suc=1; 
for ktest=1:1:length(Tests)
    time_A = cputime;
    Chosen = Tests{ktest};
    [P_Var,D,t] = Reading_Data_Base(Chosen,DB_path,Df_S);
    [ID] = Compute_slab_characteristics(P_Var.eta0DS,Df_S,n,P_Var.L0,P_Var.s0,D0,P_Var.eta0DM,Df_UM,nlm);
    % Example
    if P_Var.failed == 0
        f = 0;
        fun      = @(x) optimizing_fit(ID,nlm,x,D,t);
        f    = fminsearch(fun,f);

        ID.ID_A.n = ID.ID_A.n-abs(f);
        %ID.ID_A.fetch = 290;
        [TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
        [D,D_0D2D] = clean_data(TestData,D,t);
        res = goodnessOfFit(D,D_0D2D,'NMSE');

        
        figure(1)
        name_picture = (['AT_',num2str(ktest),'D_t']);
        plot(n.*t,D,'Color','k',LineStyle='-.',LineWidth=1.2)
        hold on
        plot(n.*t,D_0D2D,'Color','b',LineStyle='-.',LineWidth=1.2)
        title_plot = (['Test',num2str(ktest), ' $\Lambda$ = ', num2str(ID.ID_A.Lambda)]);
        title(title_plot,Interpreter="latex")
        xlim([0,10])
        ylim([0.1,1.0])
        legend('2D numerical experiment','0D numerical experiment')
        grid on
        xlabel('$\frac{n\cdot t}{t_c}$',Interpreter='latex')
        ylabel('$\frac{D}{D_0}$',Interpreter='latex')
        pt=fullfile(pt_save,name_picture);
        print(pt,'-dpng')
        clf; 

        figure(2)
        name_picture = (['BT_',num2str(ktest),'DIFF']);
        plot(n.*t,D-D_0D2D)
        xlim([0,10])
        title_plot = (['Test $D_{2D}-D_{0D}$',num2str(ktest), '$\Lambda$ = ', num2str(ID.ID_A.Lambda)]);
        title(title_plot,Interpreter="latex")

        grid on
        xlabel('$\frac{n\cdot t}{t_c}$',Interpreter='latex')
        ylabel('$\frac{D_{2D}-D_{0D}}{D_0}$',Interpreter='latex')
        pt=fullfile(pt,name_picture);
        print(pt_save,'-dpng')
        clf;

        disp(['The test @ node ', Chosen, ' show a ',num2str(res*100), ' % fit'])
        Lambda(suc)=ID.ID_A.Lambda;
        fitting_p(suc)=res;
        fetch(suc)=f; 
        L0(suc)= P_Var.L0;
        s0(suc)=P_Var.s0;
        suc = suc+1; 

    else
        disp('=====================================================================')
        disp(['You cannot optimize a failure, but, i can give you motivational tips:'])
        disp([ ' it is not always good being consistent in failing'])
        disp('So long and thanks for all the fish')
        disp('=====================================================================')
        f = nan;
        t = [];
        D = [];
        D_0D2D = [];
    end
    TB.(strcat('T',num2str(ktest))).ID = ID.ID_A;
    TB.(strcat('T',num2str(ktest))).P_Var = P_Var;
    TB.(strcat('T',num2str(ktest))).D = D;
    TB.(strcat('T',num2str(ktest))).t = t;
    TB.(strcat('T',num2str(ktest))).D0D2 = D_0D2D;
    TB.(strcat('T',num2str(ktest))).f = f;
    TB.(strcat('T',num2str(ktest))).res = res;

    B = cputime; 
    time_ktest = (B-time_A)/60;

    disp(['The test @ node ', Chosen, ' took ', num2str(time_ktest),' minutes'])

end


    figure(1)
    title_plot = ('$log_{10}(\Lambda)$ vs $fit$')
    scatter((Lambda),fitting_p.*100,"k",'filled')
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$F \%$',Interpreter='latex')
    set(gca, 'XScale', 'log')
    grid on
    pt=fullfile(pt,'SCATTER_FITTING');
    print(pt_save,'-dpng')

    figure(2)
    title_plot = ('$log_{10}(\Lambda)$ vs $n_{eff}$');
    scatter((Lambda),n-fetch,"k",'filled')
    xlabel('$log_{10}(\Lambda)$',Interpreter='latex')
    ylabel('$n_{eff}$',Interpreter='latex')
    grid on
    set(gca, 'XScale', 'log')
    pt=fullfile(pt,'SCATTER_effectiveEXP');
    print(pt_save,'-dpng')

    figure(3)
    scatter((L0/D0),n-fetch,"k",'filled')
    xlabel('$(\frac{L_0}{D_0})$',Interpreter='latex')
    ylabel('$n_{eff}$',Interpreter='latex')
    grid on
    %set(gca, 'XScale', 'log')
    pt=fullfile(pt,'SCATTERL0_effectiveEXP');
    print(pt_save,'-dpng')
    figure(4)
    scatter((L0/D0),fitting_p,"k",'filled')
    xlabel('$(\frac{L_0}{D_0})$',Interpreter='latex')
    ylabel('$F /%$',Interpreter='latex')
    grid on
    %set(gca, 'XScale', 'log')
    pt=fullfile(pt,'SCATTERL0_FITP');
    print(pt_save,'-dpng')


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

P_Var.failed = h5read(DB_path,status);

path_D = strcat(Test_name,'/TimeEvolution/Slab1D/D');
path_t = strcat(Test_name,'/TimeEvolution/time');
path_L0 = strcat(Test_name,'/Initial_Data/Initial_Condition/L0');
path_tc = strcat(Test_name,'/Initial_Data/Initial_Condition/tc');
path_s0 = strcat(Test_name,'/Initial_Data/Initial_Condition/tau0');
path_eta0DS = strcat(Test_name,'/Initial_Data/Initial_Condition/eta0');
path_eta0DM = strcat(Test_name,'/Initial_Data/Initial_Condition/Astenosphere/eta');

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
else

    % Relevant Information
    D = h5read(DB_path,path_D);
    D=D./D(1);
    t= h5read(DB_path,path_t);
    t = t./P_Var.tc;
end

end

function [res] = optimizing_fit(ID,nlm,f,D,t)
%=========================================================================$
%
%
%==========================================================================

ID.ID_A.n = ID.ID_A.n-abs(f);
%ID.ID_A.fetch = 290;
[TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
[D,D_0D2D] = clean_data(TestData,D,t);
res = goodnessOfFit(D,D_0D2D,'NMSE');
end

function [D,D_0D2D] = clean_data(TestData,D,t)

time_0D = TestData.time;
D_0D    = TestData.D_norm;
D_0D2D = interp1(time_0D,D_0D,t,'linear','extrap');
D_0D2D(isnan(D_0D2D) | D_0D2D<0.1)=0.0;
D(isnan(D) | D<0.1)=0.0;
end