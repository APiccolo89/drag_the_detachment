function [TB,FIT] = perform_optimization_DataBase(Tests,n,D0,Df_S,nlm,Df_UM,DB_path,pt_save,Res)

suc = 1;
for ktest=1:1:length(Tests)
    time_A = cputime;
    Chosen = Tests{ktest};
    [P_Var,D,t,Tau,d,tdet,topo,tauL] = Reading_Data_Base(Chosen,DB_path,Df_S,nlm);
    [ID] = Compute_slab_characteristics(P_Var.eta0DS,Df_S,n,P_Var.L0,P_Var.s0,D0,P_Var.eta0DM,P_Var.xiUM,nlm);
    % Example
    if P_Var.failed == 0
        type =1;
        ID.ID_A.cut_off_Mantle = 1.0;
        ID.ID_A.cut_off_Slab   = 1.0;
        ID.ID_A.iteration      = 0;
        ID.ID_A.flag =0;
        if ktest == 1
            test_Gif = 1;
        else
            test_Gif = 0;
        end
        f = [abs(double(-d/2)),1];
        fun      = @(x) optimizing_fit(ID,nlm,x,D,t,tdet,type,test_Gif);
        if type == 2
            f    = fsolve(fun,f);
            ID.ID_A.flag = 1;
            ID.ID_A.fetch(1) = abs(f(1));
            ID.ID_A.fetch(2) = abs(f(2));
            [TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
            [D,D_0D2D,t] = clean_data(TestData,real(D),t);
            res = (tdet-TestData.t_det.*ID.n)./tdet;
        else
            f    = fminsearch(fun,f);
            ID.ID_A.flag = 1;
            ID.ID_A.fetch(1) = abs(f(1));
            ID.ID_A.fetch(2) = abs(f(2));
            [TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
            [D,D_0D2D,t] = clean_data(TestData,real(D),t);
            res = goodnessOfFit(D,D_0D2D,'NMSE');

        end

        disp(['The test @ node ', Chosen, ' show a ',num2str(res*100), ' % fit'])

        if nlm.islinear == 0
            FIT.Lambda(suc)=ID.ID_A.Lambda./P_Var.xiUM;
        else
            FIT.Lambda(suc)=ID.ID_A.Lambda;
        end
        FIT.fitting_p(suc)=res;
        FIT.fetch(1,suc)=f(1);
        FIT.fetch(2,suc)=f(2);
        FIT.eta0DM(suc)=P_Var.eta0DM;
        FIT.L0(suc)= P_Var.L0;
        FIT.s0(suc)=P_Var.s0;
        FIT.Dp(suc) = d;
        FIT.Detachment(1,suc) = n*TestData.t_det;
        FIT.Detachment(2,suc) = tdet;
        FIT.Res(suc)= Res(ktest);
       
        tau0D = TestData.tau;
        t0D   = TestData.time;
       
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

        tau0D =[];
        t0D   = [];
       
    end
    TB.(strcat('T',num2str(ktest))).ID = ID.ID_A;
    TB.(strcat('T',num2str(ktest))).tc = ID.tc;

    TB.(strcat('T',num2str(ktest))).P_Var = P_Var;
    TB.(strcat('T',num2str(ktest))).D = D;
    TB.(strcat('T',num2str(ktest))).t = t;
    TB.(strcat('T',num2str(ktest))).t0D =t0D;
    TB.(strcat('T',num2str(ktest))).tau0D= tau0D;
    TB.(strcat('T',num2str(ktest))).tau = Tau; 
    TB.(strcat('T',num2str(ktest))).topo = topo;
    TB.(strcat('T',num2str(ktest))).D0D2 = D_0D2D;
    TB.(strcat('T',num2str(ktest))).f = f;
    TB.(strcat('T',num2str(ktest))).res = res;

    B = cputime;
    time_ktest = (B-time_A)/60;
        t = [];
        D_0D2D = [];
        D_0D2DW = [];
        Tau02D = [];
        error = [];
    disp(['The test @ node ', Chosen, ' took ', num2str(time_ktest),' minutes'])

end

%close all;
%clf;
end