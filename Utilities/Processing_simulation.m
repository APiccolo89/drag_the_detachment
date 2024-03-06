function [Temp]=Processing_simulation(eta0DS,Df_S,n,l0,s0,D0,eta0DM,Benchmark,Df_UM,nlm)
    % Produce the structure containing all the relevant data for running
    % the simulation
    [ID] = Compute_slab_characteristics(eta0DS,Df_S,n,l0,s0,D0,eta0DM,Df_UM,nlm);
    disp(['log10(Psi) = ',num2str(log10(ID.Psi))]); 
    disp('Test featuring log10 Psi_c >0 are not performed, and are considered failed')
    if log10(ID.Psi/(1+ID.Df_UM.*ID.ID_A.tau_mc.^(ID.n-1)))>0 && nlm.Linear == 0 
        Testdata = [];
        Testdata.Failed = 2; % (0! Success, 1! Failed, 2! not worth of being even born)
    else
        Testdata = Run_Simulation_DragA(ID.ID_A,nlm);
    end
    %Store the data, such that is possible to understand which portion of
    %the parametric space is critical. 
    Temp = Testdata;
    Temp=Testdata;
    Temp.initial_data = ID;
    ID = [];
    Testdata = [];
    % Empty the trashbin.
end

