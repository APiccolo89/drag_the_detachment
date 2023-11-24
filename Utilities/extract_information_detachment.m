function [Data_S] = extract_information_detachment(S,hyp,nlm)

if hyp == 1 %check if it is a structure that has higher hierarchy
    fn = fieldnames(S);
    % Prepare Data_S array
    L_test = 0;
    for k = 1:numel(fn)
        fn_l=(fieldnames(S.(fn{k})));
        L_test = L_test+length(fn_l);
    end

else
    fn = [];
    fn_l = fieldnames(S);
    L_test = length(fn_l);
end


%function to setup structure:
Data_S = Detachment_DB;
Data_S.Create_Vectors(L_test);

i=1;
if isempty(fn)
    [Data_S] = extract_raw_data(Data_S,S,fn_l,i,nlm);
else

    for k = 1:numel(fn)
        S_1=S.(fn{k});
        [Data_S,i] = extract_raw_data(Data_S,S_1,fn_l,i,nlm);
    end
end

end

function [Data_S,i]=extract_raw_data(Data_S,S,fn_l,i,nlm)

F = fieldnames(S);
for k = 1:numel(F)

    TD = S.(F{k});
    Data_S.Failed(i) = TD.Failed; 
    if TD.Failed ~= 2
        Data_S.max_step(i) = TD.option.max_step;
        Data_S.iteration(i) = TD.option.iteration;
    end
    Data_S.Lambda(i) = TD.initial_data.Lambda;
    Data_S.Psi(i) = (TD.initial_data.Psi);
    Data_S.xiUS(i)= TD.initial_data.Df_S;
    Data_S.L0(i) = TD.initial_data.l0;
    Data_S.tc(i) = TD.initial_data.tc;
    Data_S.n(i) = TD.initial_data.n;
    Data_S.eps_c(i) = 1./TD.initial_data.tc;
    Data_S.eta0DS(i) = TD.initial_data.eta0DS;
    Data_S.eta0DUM(i) = TD.initial_data.eta0DM;
    Data_S.tau0(i)   = TD.initial_data.s0;
    Data_S.BnUM(i)   = TD.initial_data.B_n_um;
    Data_S.BdUM(i)   = TD.initial_data.B_d_um; 
    Data_S.BnS(i)   = TD.initial_data.B_n;
    Data_S.BdS(i)   = TD.initial_data.B_d; 

    if nlm.islinear ==0
        Data_S.xiUM(i) = TD.initial_data.Df_UM;
    end
    if isfield(TD,'t_det')==1
    if ~isempty(TD.t_det) && isreal(TD.D_norm)
        Data_S.tdet(i) = TD.initial_data.n*TD.t_det;
        Data_S.tau_max(i)= TD.t_t_max;
        Data_S.time_tau_max(i) = TD.time_t_M(end)*TD.initial_data.n;
        Data_S.tau_det(i)= TD.t_t_det;
        Data_S.tau_real_initial(i) = TD.tau(3,1); 
        Data_S.tau_drag_initial(i) = -TD.tau(2,1);
        % compute the tc assuming that the reference stress is the initial
        % one: 
        tau_eff_in = Data_S.tau_real_initial(i)*TD.initial_data.s0;
        eps_c_drag = TD.initial_data.B_d.*tau_eff_in+TD.initial_data.B_n.*tau_eff_in.^TD.initial_data.n;
        Data_S.tc_drag(i) = 1./eps_c_drag; 
        Data_S.Lambda0(i) = TD.Lambda(1);

        if nlm.iteration > 0
            Data_S.xdisl(i)   = TD.eps_mantle(1,1)./TD.eps_mantle(3,1); 
            Data_S.tau_um_0(i) = TD.tau_mantle(1);
            Data_S.eps_um_0(i) = TD.eps_mantle(3,1);
            Data_S.eta_um_0(i) = TD.eta_um(1);
            Data_S.tau_mc(i) = TD.initial_data.ID_A.tau_mc;
        end
       Data_S.dDdt_0(i) = TD.dDdt(1); 
       
    else
        Data_S.tdet(i) = nan;
        Data_S.tau_max(i)= nan;
        Data_S.time_tau_max(i) = nan;
        Data_S.tau_det(i)= nan;
        Data_S.tau_real_initial(i) = nan; 
        Data_S.tc_drag(i) = nan; 
        Data_S.tau_drag_initial(i) = nan; 
        Data_S.Lambda0(i)= nan;
        if nlm.iteration > 0
            Data_S.tau_um_0(i) = nan;
            Data_S.xdisl(i)   = nan; 
            Data_S.eps_um_0(i) = nan;
            Data_S.eta_um_0(i) = nan;
            Data_S.tau_mc(i) = nan;

        end

    end
    if isfield(TD,'Meta_dataReal')
        MD = TD.Meta_dataReal; 
        fname = fieldnames(MD);
        for in = 1:numel(fname)
            Data_S.(fname{in})(i)=MD.(fname{in});
        end

    end
    else
        Data_S.tdet(i) = nan;
        Data_S.tau_max(i)= nan;
        Data_S.time_tau_max(i) = nan;
        Data_S.tau_det(i)= nan;
        Data_S.tau_real_initial(i) = nan; 
        Data_S.tc_drag(i) = nan; 
        Data_S.tau_drag_initial(i) = nan; 
    end
        i = 1+i;
end
end

