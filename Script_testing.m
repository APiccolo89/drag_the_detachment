%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter to test: 
% To do: 
% Generate a data structure and labelling. 

l0_v    = (300e3);                 % initial length
D0      = 80e3;                  % thickness

% Slab Rheology
eta0_v = [1e21];                   % [Pas] refernce power law viscosity slab 
Df   = 10;                     % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress 
n     = 3.5;                   % [n.d.] pre-exponential factor
s0_v  = (100e6);                 % [Pa]  reference buoyancy stress

% Upper Mantle
etaum_v =10.^(16:0.1:20);    % [Pa.s]vector of the mantle viscosity



% Function to run the ensamble of test 
[Testdata] = Run_Simulations(D0,l0_v,eta0_v,Df,n,s0_v,etaum_v);


% Function to plot the results 



%========================================================================%
% Function utilities                                                     %
%========================================================================%
function [Testdata] = Run_Simulations(D0,l0_v,eta0_v,Df,n,s0_v,etaum_v)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Output:
     % 
     % [Testdata] => structure containing the results of each of the
     %               simulation performed 
     % 
     % Input: 
     % 
     % D0,l0_v => Initial thickness and length (l0_v vector of data)
     %
     % eta0_v => reference power law viscosity at reference stress (tau0)
     %       n power law exponent [vector]
     % 
     % Df    =>viscosity contrast between diffusion/dislocation creep @
     %       reference stress (tau0)
     %
     % etaum_v => Viscosity of the upper mantle [vector]. 
     %
     % Alg:
     % -> loop over the data and generate combination of simulation 
     % -> construct the characteristic values 
     % -> run simulation per each combination of parameter 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     ium   = length(etaum_v);    % number of upper mantle viscosities
     iut0  = length(s0_v) ;      % number of reference buoyancy stresses
     iuet0 = length(eta0_v);     % number of reference viscosities of the slab
     iuL0  = length(l0_v)   ;    % number of tested initial length of the slab

     if length(n)>1 || length(D0)>1 
        error('You cannot test more than one power law exponent or initial thickness at once')
 
     end

  %loop to construct the initial data to run a simulation (to change into a
  %more performant function)
  for i=1:ium
      etaum = etaum_v(i);
      for j = 1:iut0
          s0 = s0_v(j);
        for k = 1:iuet0
            eta0 = eta0_v(k);
           for z = 1:iuL0
                l0 = l0_v(z);
                % Compute the characteristic for the simulation 
                [drho,B_d,B_n,tc,ec] = Compute_slab_characteristics(eta0,Df,n,l0,s0);
                % Create a labeling for the simulation
                str0   = strcat('Sim_D0_',num2str(int16(D0/1e3)));
                str1   = strcat('em_',num2str(i));
                str2   = strcat('s0_',num2str(j));
                str3   = strcat('eta0_',num2str(k));
                str4   = strcat('L0_',num2str(z));
                l_simulation = strcat(str0,str1,str2,str3,str4);
                disp(['Simulation ',l_simulation, 'is starting'])
                % Save into a simple data structure the initial data of
                % the simulation
                string_ID = {'B_d','B_n','s0','n','eta0','Df','drho','D0','l0','etaum','tc','ec'};
                
                for is = 1:numel(string_ID)
                    ID.(string_ID{is}) = eval(string_ID{is});
                end
                % Run a simulation with a specific combination of parameter               
                Testdata = Run_Simulation_Drag(ID);
                % 
                Tests.(l_simulation) = Testdata;
                Tests.(l_simulation).initial_data = ID; 
                ID = [];
                Testdata = []; 
 
           end
        end
      end
  end

end







function [drho,B_d,B_n,tc,ec] = Compute_slab_characteristics(eta0,Df,n,l0,s0)
        % Output: 
        % drho = density contrast
        % B_d  = diffusion creep pre-exponential factor
        % B_n  = dislocation creep pre-exponential factor
        % tc   = characteristic time (s)
        % ec   = characterstistc strain rate (1/s)
        % Input
        % eta0 = reference viscosity at referenc stress
        % s0   = reference stress
        % l0   = initial lenght
        % D0   = initial thickness
        % n    = power law exponent
        % Df   = viscosity contrast between diffusion and dislocation at reference
        %        stress
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        drho   = 2*s0/(9.81*l0);     % delta rho 
        B_n   = s0^(1-n)/eta0;       % compliance dislocation
        B_d   = 1/(Df*eta0);         % compliance diffusion
        ec    = (B_n*s0^n+B_d*s0);   % characteristic strain rate
        tc    = 1/ec             ;   % characteristic time scale 

end




% 
% 
% 
% for i = 1:length(etaUM_vec)
%     etaUM = etaUM_vec(i); 
%     Testdata = Compute_dddt_Drag(D,B_n,B_d,n,L0,D0,drho,etaUM,tau0,name_1);
%     name_2 = append('T',name_1,num2str(i));
%     Tests.(name_2) = Testdata;
%     Testdata = [];
%     eta_UM = [];  
% end
% 
% 
% 
% 
% % Plot function
% i = 1 ; 
% fn = fieldnames(Tests);
% cc = jet(length(etaUM_vec));
% for k = 1:numel(fn) 
% 
% TD = Tests.(fn{k});
% T = TD(1,:);
% D = TD(2,:);
% hold on 
% plot(T,D,'Color',cc(i,:))
% 
% %set(gca, 'YScale', 'log')
% grid on 
% xlim([0,10])
% ylim([10^(-2),10^(0)])
% 
% xlabel('t/tc [n.d.]')
% 
% i = i+1; 
% end
% print('Global_Test','-dpng')
% clf; 
% close; 