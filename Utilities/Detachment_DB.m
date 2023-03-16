classdef Detachment_DB < handle
    %This classes is a datatype structure that allows to store data of the
    %simulation that are relevant and creates the automatically. 
    % Methods 1: initializate the vectors (input the number of tests)

    properties
        L0           % Initial Length
        D0           % Initial Thickness
        tau0         % Reference Stress
        n            % Stress exponent
        eta0DS       % Reference diffusion creep viscosity Slab
        eta0DUM      % Reference diffusion creep viscosity Upper Mantle
        xiUM         % Viscosity contrast between diffusion and dislocation at reference condition of the upper mantle
        xiUS         % Viscosity contrast between diffusion and dislocation at reference condition of the Slab
        eps_c        % Characteristic strain rate
        Lambda       % Initial Lambda
        Psi          % Viscosity contrast at reference condition between mantle and slab 
        s            % Length scale of the convection
        tc           % Characteristic timescale 
        tc_drag      % Characteristic timescale using the drag stress
        tdet         % Detachment time (t/tc)
        tau_det      % tau at the detachment 
        tau_max      % tau max 
        time_tau_max % time at which the tau max is happening
        tau_real_initial % initial stress at the initial stage
        tau_drag_initial % initial drag stress 
        T_Slab % Field that is useful for the realistic case: the average temperature of the slab in Kelvin
        Tp    % Mantle potential temperature
        Vn    % Activaction volume 
        Vd    % Activaction volume 
        Cd    % Exponential factor diffusion
        Cn    % Exponential factor dislocation
        phi   % adiabatic gradient
        w     % weight force
        Pr    % Reference pressure
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