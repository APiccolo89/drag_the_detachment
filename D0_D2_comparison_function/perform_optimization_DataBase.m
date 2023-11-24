function [TB,FIT] = perform_optimization_DataBase(Tests,n,D0,Df_S,nlm,Df_UM,DB_path,pt_save,Res,number_fetch,matlab_version)

suc = 1;

for ktest=1:1:length(Tests)

    time_A = cputime;

    Chosen = Tests{ktest};

    [P_Var,D,t,Tau,d,tdet,topo,tauL,d2,D_matrix] = Reading_Data_Base(Chosen,DB_path,Df_S,nlm);

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
        if number_fetch == 1

            f = [abs(double(d./2))];
        else

            f = [abs(double(d./2)),abs(1)];
        end

        fun      = @(x) optimizing_fit(ID,nlm,x,D,t,tdet,type,test_Gif,number_fetch,matlab_version);
        res_gf = [];
        if type == 2
            f    = fsolve(fun,f);
            ID.ID_A.flag = 1;
            for i = 1:number_fetch
                ID.ID_A.fetch(i) = abs(f(i));
            end
            [TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
            [D,D_0D2D,t] = clean_data(TestData,real(D),t);
            res_gf = goodnessOfFit(D,D_0D2D,'NMSE');

            res = (tdet-TestData.t_det.*ID.n)./tdet;

        else
            f    = fminsearch(fun,f);
            ID.ID_A.flag = 1;
            for i = 1:number_fetch
                ID.ID_A.fetch(i) = abs(f(i));
            end
            [TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
            [D,D_0D2D,t] = clean_data(TestData,real(D),t);
            res = goodnessOfFit(D,D_0D2D,'NRMSE');
            res_gf = res;
        end
        if matlab_version<=2020
            disp(['The test @ node ', Chosen, ' show a ',num2str((1-res_gf)*100), ' % fit'])
        else
            disp(['The test @ node ', Chosen, ' show a ',num2str((res_gf)*100), ' % fit'])

        end

        if nlm.islinear == 0
            FIT.Lambda(suc)=ID.ID_A.Lambda./(1+P_Var.xiUM.*ID.ID_A.tau_mc.^(ID.n-1));
        else
            FIT.Lambda(suc)=ID.ID_A.Lambda;
        end
        FIT.fitting_p(suc)=res;
        FIT.fetch(1,suc)=f(1);

        if number_fetch==1
            FIT.fetch(2,suc)=1.0;

        else
            FIT.fetch(2,suc)=f(2);

        end
        FIT.eta0DM(suc)=P_Var.eta0DM;
        FIT.xium(suc) = P_Var.xiUM;
        FIT.xius(suc) = 10.0;
        FIT.eta0DS(suc)  = P_Var.eta0DS;
        FIT.L0(suc)= P_Var.L0;
        FIT.s0(suc)=P_Var.s0;
        FIT.Dp(suc) = d;
        FIT.Dp2(suc) = d2;
        FIT.Detachment(1,suc) = n*TestData.t_det;
        FIT.Detachment(2,suc) = tdet;
        FIT.Res(suc)= Res(ktest);
        FIT.res_gf(suc)=res_gf;
        Detachment_0D = n*TestData.t_det;
        Detachment_2D = tdet;
        tau0D = TestData.tau;
        t0D   = TestData.time;

        suc = suc+1;
  
        if number_fetch >1
            figure(3);
            ax = gca;
            ax.XColor = [0,0,0];
            ax.YColor = [0,0,0];
            ax.LineWidth = 1.2;
            ax.Box     = 'on';
            filename=(['T_2f',num2str(ktest),'.png']);
            pt = fullfile(pt_save,filename);
            print(pt,'-dpng','-r0')
            clf;
        end
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
        Detachment_0D = [];
        Detachment_2D = [];

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
    TB.(strcat('T',num2str(ktest))).D_matrix = D_matrix;
    TB.(strcat('T',num2str(ktest))).f = f;
    TB.(strcat('T',num2str(ktest))).res = res;
    TB.(strcat('T',num2str(ktest))).Detachment_0D = Detachment_0D;
    TB.(strcat('T',num2str(ktest))).Detachment_2D = Detachment_2D;
    TB.(strcat('T',num2str(ktest))).tauL = tauL;
    TB.(strcat('T',num2str(ktest))).Dp2 = d2;
    if isempty(f)
        TB.(strcat('T',num2str(ktest))).Dp = nan;

    else
        TB.(strcat('T',num2str(ktest))).Dp = -f(1)*2*P_Var.L0./1e3;
    end
    B = cputime;
    time_ktest = (B-time_A)/60;
    t = [];
    D_0D2D = [];
    D_0D2DW = [];
    Tau02D = [];
    error = [];
    disp(['The test @ node ', Chosen, ' took ', num2str(time_ktest),' minutes'])

end
end





