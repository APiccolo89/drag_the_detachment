classdef Detachment_DB < handle
    %This classes is a datatype structure that allows to store data of the
    %simulation that are relevant and creates the automatically. 
    % Methods 1: initializate the vectors (input the number of tests)

    properties
        L0           % Initial Length[Input]
        D0           % Initial Thickness[Input]
        tau0         % Reference Stress[Input]
        n            % Stress exponent[Input]
        eta0DS       % Reference diffusion creep viscosity Slab[Input]
        eta0DUM      % Reference diffusion creep viscosity Upper Mantle[Input]
        xiUM         % Viscosity contrast between diffusion and dislocation at reference condition of the upper mantle[Input]
        xiUS         % Viscosity contrast between diffusion and dislocation at reference condition of the Slab[Input]
        eps_c        % Characteristic strain rate[Input]
        Lambda       % Initial Lambda[Input]
        Psi          % Viscosity contrast at reference condition between mantle and slab [Input]
        s            % Length scale of the convection[Input]
        tc           % Characteristic timescale [Input]
        tc_drag      % Characteristic timescale using the drag stress[numeric]
        tdet         % Detachment time (t/tc)[numeric]
        tau_det      % tau at the detachment[numeric] 
        tau_max      % tau max [numeric]
        time_tau_max % time at which the tau max is happening[numeric]
        tau_real_initial % initial stress at the initial stage[numeric]
        tau_drag_initial % initial drag stress [numeric]
        T_Slab % Field that is useful for the realistic case: the average temperature of the slab in Kelvin[Input]
        Tp    % Mantle potential temperature[Input]
        Vn    % Activaction volume [Input]
        Vd    % Activaction volume [Input]
        Cd    % Exponential factor diffusion[Input]
        Cn    % Exponential factor dislocation[Input]
        phi   % adiabatic gradient[Input]
        w     % weight force[Input]
        Pr    % Reference pressure[Input]
        Lambda0 % initial Lambda [Numeric]
        BdUM % Diffusion creep exponential factor[Input]
        BnUM % Dislocation creep exponential factor[Input]
        BdS % Diffusion creep exponential factor[Input]
        BnS % Dislocation creep exponential factor [Input]
        tau_um_0 % initial stress of the upper mantle [Numeric]
        xdisl % initial partition of dislocation diffusion creep [Numeric]
        xdisl_mean %mean value of xdisl 
        xdisl_max %max value of xdisl 
        eps_um_0 %initial strain rate of the upper mantle [Numeric]
        dDdt_0 %initial dDdt at the beginning of the simulation [Numeric]
        eta_um_0 % initial upper mantle viscosity [Numeric]
        tau_mc % characteristic mantle length scale [Numeric/Input]
        max_step % resolution of the timestep implicit model
        iteration % how many time do I need to use the quick hack
        Failed %flag of the test
    end

    methods
        function obj = Create_Vectors(obj,n_tests)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            p = properties(obj);
            for iprop = 1:length(p)
                obj.(p{iprop}) = zeros(n_tests,1);
            end
        end
        
    end
end