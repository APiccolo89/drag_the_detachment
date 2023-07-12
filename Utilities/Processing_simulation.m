function [Temp]=Processing_simulation(eta0DS,Df_S,n,l0,s0,D0,eta0DM,Benchmark,Df_UM,nlm)
    % Produce the structure containing all the relevant data for running
    % the simulation
    [ID] = Compute_slab_characteristics(eta0DS,Df_S,n,l0,s0,D0,eta0DM,Df_UM,nlm);
    %ID.ID_A.cut_off_Mantle = 1.0;
    %ID.ID_A.cut_off_Slab   = 1.0;    % Run a simulation with a specific combination of parameter
    %Testdata = Run_Simulation_Drag(ID,Benchmark,nlm);
    disp(['log10(Lambda) = ',num2str(log10(ID.ID_A.Lambda))]); 
    Testdata = Run_Simulation_DragA(ID.ID_A,nlm);
    %interpolation routine from AD->D
    Temp = Testdata;
    Temp=Testdata;
    Temp.initial_data = ID;
    ID = [];
    Testdata = [];
end

