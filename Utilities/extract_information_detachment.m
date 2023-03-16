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

    Data_S.Lambda(i) = TD.initial_data.Lambda;
    Data_S.Psi(i) = (TD.initial_data.Psi);
    Data_S.xiUS(i)= TD.initial_data.Df_S;
    Data_S.L0(i) = TD.initial_data.l0;
    Data_S.tc(i) = TD.initial_data.tc;
    Data_S.n(i) = TD.initial_data.n;
    Data_S.eps_c(i) = 1./TD.initial_data.tc;
    Data_S.eta0DS(i) = TD.initial_data.eta0DS;
    Data_S.eta0DUM(i) = TD.initial_data.eta0DM;
    if nlm.islinear ==0
        Data_S.xiUM(i) = TD.initial_data.Df_UM;
    end

    if ~isempty(TD.t_det) && isreal(TD.D_norm)
        Data_S.tdet(i) = TD.initial_data.n*TD.t_det;
        Data_S.tau_max(i)= TD.t_t_max;
        Data_S.time_tau_max(i) = TD.time_t_M*TD.initial_data.n;
        Data_S.tau_det(i)= TD.t_t_det;
        Data_S.tau_real_initial(i) = TD.tau(3,1); 
        Data_S.tau_drag_initial(i) = -TD.tau(2,1);
        % compute the tc assuming that the reference stress is the initial
        % one: 
        tau_eff_in = Data_S.tau_real_initial(i)*TD.initial_data.s0;
        eps_c_drag = TD.initial_data.B_d.*tau_eff_in+TD.initial_data.B_n.*tau_eff_in.^TD.initial_data.n;
        Data_S.tc_drag(i) = 1./eps_c_drag; 
        

    else
        Data_S.tdet(i) = nan;
        Data_S.tau_max(i)= nan;
        Data_S.time_tau_max(i) = nan;
        Data_S.tau_det(i)= nan;
        Data_S.tau_real_initial(i) = nan; 
        Data_S.tc_drag(i) = nan; 
        Data_S.tau_drag_initial(i) = nan; 
    end
    if isfield(TD,'Meta_dataReal')
        MD = TD.Meta_dataReal; 
        fname = fieldnames(MD);
        for in = 1:numel(fname)
            Data_S.(fname{in})(i)=MD.(fname{in});
        end

    end
    i = 1+i;
end
end

