%==========================================================================
clear all
close all
clf
%==========================================================================
addpath Adimensional\
addpath Dimensional\
addpath Utilities\
Benchmark = 0.0 ; % Benchmark activaction flag
%==========================================================================
%Parameter to test: 
%==========================================================================
% Testing parameter: 
%==========================================================================
% l0_v = vector containing all the lenght to test
% D0   = initial thickness 
% s0_v = buoyancy stress applied 
%==========================================================================
% eta0DS_v = [vector] reference linear viscosity of the slab.
% eta0DM_v = [vector] reference linear viscosity of the mantle.
% DfS_v    = [vector] viscosity contrast between power-law rheology and
% linear one at reference stress condition
% DfM_v    = [vector] viscosity contrast between power-law rheology and
% linear one at reference condition [deactivated if upper mantle is linear]
% n_v      = [vector] power law exponent {same Slab-Mantle}
% nlm      = activation flag non linear mantle. 
% These are the primary values used to construct the set of experiments.
% Instead of using directly Lambda as primary value, I decided to use
% Earth-like range parameters that are used to compute the relevant data
% such that is possible to associated set of parameters to the
% dimensionless group. 
%==========================================================================
%[OUTPUT]: => [Tests.mat] Database featuring relevant data for each tests
%i.e. Table_Experiments:     T1 => [D1-------Dn|Detachment data]
%     Numerical_Experiments: T1 -
%                               |_ INITIAL DATA STRUCTURE
%                               |_ D_norm [1,nstep] 
%                               |_ Stress [3,nstep]: 1:Effective stress 2:
%                                   Bouyancy stress 3: Drag integral stress
%                               |_ Def [3, nstep]:  1: total strain; 2:
%                                   dislocation creep strain; 3: Diffusion
%                                   creep strain
%                               |_ Lambda[1,nstep]
%                               |_ UM_vis[1,nstep]
%                               |_ time_Vector[1,nstep]
%                            Tn    
% The numbering of the test depends on the layout of the matrix created in 
% f. Run Simulations[]
% Pictures: saved outside the actual folder of the script. And divided
% using the stress exponent.
%==========================================================================
% Some important note: The numerical code has not a default viscosity cut
% off to avoid the limitation that usual 2D numerical code has. On the
% other hand, there are some combination of parameter that are not working
% well together, and this is a consequence that the non linear viscosity of
% of the slab is too unconstrained. xiUS must be used carefully, due to
% these reason. 
%==========================================================================
% Folder output: 
ptsave = ['../Tests_Results_NLINEAR']; 
% if the folder does not exist, create the folder. 
if not(isfolder(ptsave))
    mkdir(ptsave)
end
nlm       = Problem_type.Linear;   % Switching the position of linear-non_linear activate the non linear upper mantle routine. 
if nlm.islinear == 1
    ptsave = ['../Tests_Results_LINEAR']; 

else
    ptsave = ['../Tests_Results_NLINEAR']; 

end
xiUM_v = [];
% Linear Test Manuscript 
if nlm.islinear == 1
    %Geometric properties of the Slab 
    l0_v      = (300e3:50e3:600e3);         % initial length [m]
    D0_v      = 80e3;                       % thickness [m]
    s0_v      = [60e6:10e6:240e6];        % reference buoyancy stress [Pa]
    % Slab Rheology
    eta0DS_v  = 10.^[22:0.5:24];                      % [Pas] reference diffusion creep viscosity of the slab 
    xiUS_v     = [10.0];               % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress 
    n_v       = [2.0,3.5,5.0];            % [n.d.] pre-exponential factor
    % Upper Mantle
    eta0DM_v  = 10.^[18:0.5:21.5];                       % [Pas] reference diffusion creep viscosity of the upper mantle 
end
if nlm.islinear==0  
    %Geometric properties of the Slab 
    l0_v      = (300e3:50e3:600e3);         % initial length [m]
    D0_v      = 80e3;                      % thickness [m]
    s0_v      = [60e6:10e6:240e6];        % reference buoyancy stress [Pa]
    % Slab Rheology
    eta0DS_v  = 10.^[22,23,24];                      % [Pas] reference diffusion creep viscosity of the slab 
    xiUS_v     = [1.0,100.0];               % [n.d.]viscosity contrast between diffusion and dislocation creep at the reference stress 
    n_v       = [2.0,3.5,5.0];            % [n.d.] pre-exponential factor
    % Upper Mantle
    eta0DM_v  = 10.^[18.0:1.0:23];                       % [Pas] reference diffusion creep viscosity of the upper mantle 
    xiUM_v     = 10.^[2,3,4,5];                  % [n.d.] viscosity contrast between diffusion and dislocation creep at reference stress (UM)
end


plot = 1; 
%==========================================================================

for i=1:length(n_v)
    n = n_v(i);
    name_tests= strcat('Tests_n_',num2str(n));
    % Function to run the ensamble of test 
    [Tests] = Run_Simulations(D0_v,l0_v,eta0DS_v,xiUS_v,n,s0_v,eta0DM_v,xiUM_v,Benchmark,nlm);
    if plot == 1 
        plot_results(Tests,name_tests,ptsave,nlm)
    end
    T.(strcat('n_',num2str(floor(n))))= Tests;    
    Tests = []; 
end
%%
[Data_S] = extract_information_detachment(T,1,nlm);

if nlm.islinear==0
    figure_nonlinearmanuscript(Data_S);
else
    figure_linear_manuscript(Data_S);
end
%%




function figure_nonlinearmanuscript(Data_S)



%%
figure(1)
s_ = Data_S.xiUS;
x = Data_S.Lambda./(1+Data_S.xiUM);
y = ((Data_S.tdet));
y = 1./y; 
z = Data_S.eta0DS; 
scatter(x(s_==1.0 & Data_S.n==3.5),y(s_==1.0 & Data_S.n==3.5),10,'d','filled','b','MarkerEdgeColor','k');
hold on 
scatter(x(s_==100.0 & Data_S.n==3.5),y(s_==100.0 &Data_S.n==3.5),10,'square','filled','r','MarkerEdgeColor','k');
legend({'$ \xi^{S} = 1.0$','$ \xi^{S} = 100.0 $'},'Interpreter','latex');
grid on
box on
set(gca, 'XScale', 'log')
grid on
box on
xlim([10^(-8) 10^0])
ylim([-0.1 1.1])
print('../FigureNLM1','-dpng')


figure(2)
s_ = Data_S.xiUS;
x = Data_S.Lambda./(1+Data_S.xiUM);
y = ((Data_S.tau_real_initial));
z = Data_S.eta0DS;
scatter(x(s_==1.0 & Data_S.n==3.5),y(s_==1.0 & Data_S.n==3.5),10,'d','filled','b','MarkerEdgeColor','k');
hold on 
scatter(x(s_==100.0 & Data_S.n==3.5),y(s_==100.0 &Data_S.n==3.5),10,'square','filled','r','MarkerEdgeColor','k');
legend({'$ \xi^{S} = 1.0$','$ \xi^{S} = 100.0 $'},'Interpreter','latex'); 
grid on
box on
set(gca, 'XScale', 'log')
grid on
box on
set(gca, 'XScale', 'log')
xlim([10^(-8) 10^0])
 ylim([-0.1,1.1])
print('../FigureNLM2','-dpng')

figure(3)
s_ = Data_S.xiUS;
x = Data_S.Lambda./(1+Data_S.xiUM);
y = ((Data_S.tau_det));
z = Data_S.eta0DS;
scatter(x(s_==1.0 & Data_S.n==3.5),y(s_==1.0 & Data_S.n==3.5),10,'d','filled','b','MarkerEdgeColor','k');
hold on 
scatter(x(s_==100.0 & Data_S.n==3.5),y(s_==100.0 &Data_S.n==3.5),10,'square','filled','r','MarkerEdgeColor','k');
legend({'$ \xi^{S} = 1.0$','$ \xi^{S} = 100.0 $'},'Interpreter','latex');
grid on
box on
set(gca, 'XScale', 'log')
xlim([10^(-8) 10^0])
print('../FigureNLM3','-dpng')

end
function figure_linear_manuscript(Data_S)
% Figure S3 Manuscript 
n_p  = 3.5; 
n    = Data_S.n;
x = Data_S.Lambda(n==n_p);
y = (Data_S.tdet(n==n_p));
z = Data_S.xiUS(n==n_p);
y = 1./y; 
% plot effects xi, for a given n. 
% Select data
i_xi  = unique(z); 
markers_s = ["square","o","^","x"];
colors_s  = ["#7E2F8E","#EDB120","#A2142F","#77AC30"];
figure(1)
hold on 
for i = 1:length(i_xi)
    s=scatter(x(z==i_xi(i)),y(z==i_xi(i)),20);
    s.MarkerEdgeColor = 'k';
    s.Marker          = markers_s{i};
    s.MarkerFaceColor = colors_s{i};
    s.LineWidth       = 0.2; 
end
legend({'$ \xi^{S} = 1.0$','$ \xi^{S} = 10.0 $','$ \xi^{S} = 100.0 $','$\xi^{S} = 1000.0$'},'Interpreter','latex');
grid on
box on
set(gca, 'XScale', 'log')
xlim([10^(-8),5])
ylim([-0.1,1.1])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
print('../FigureS3_A','-dpng')
clf;

x = [];
y = [];
z = []; 

xi  = 10.0; 
x = Data_S.Lambda(Data_S.xiUS==xi);
y = Data_S.n(Data_S.xiUS==xi)./Data_S.tdet(Data_S.xiUS==xi);
z = Data_S.n(Data_S.xiUS==xi);
i_xi  = unique(z); 

figure(2)
hold on 
for i = 1:length(i_xi)
    s=scatter(x(z==i_xi(i)),y(z==i_xi(i)),20);
    s.MarkerEdgeColor = 'k';
    s.Marker          = markers_s{i};
    s.MarkerFaceColor = colors_s{i};
    s.LineWidth       = 0.2; 
end
legend({(['$ n  =$',num2str(i_xi(1))]),(['$ n  =$',num2str(i_xi(2))]),(['$ n  =$',num2str(i_xi(3))])},'Interpreter','latex');
grid on
box on
xlim([10^(-8),5])
ylim([-0.1,7.5])
set(gca, 'XScale', 'log')

%set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
print('../FigureS3_B','-dpng')
clf;

figure(3)
y = Data_S.n(Data_S.xiUS==xi).*Data_S.time_tau_max(Data_S.xiUS==xi);

hold on 
for i = 1:length(i_xi)
    s=scatter(x(z==i_xi(i)),y(z==i_xi(i)),20);
    s.MarkerEdgeColor = 'k';
    s.Marker          = 'o';
    s.MarkerFaceColor = colors_s{i};
    s.LineWidth       = 0.2; 
end
legend({(['$ n  =$',num2str(i_xi(1))]),(['$ n  =$',num2str(i_xi(2))]),(['$ n  =$',num2str(i_xi(3))])},'Interpreter','latex');
grid on
box on
xlim([10^(-8),5])
ylim([-0.1,10.0])
set(gca, 'XScale', 'log')
%set(gca,'XTickLabel',[])
%set(gca,'YTickLabel',[])
print('../FigureS3_C','-dpng')

end
%%
% % Figure S3 Manuscript 
% n_p  = 3.5; 
% x = Data_S(1,Data_S(7,:)==n_p);
% y = 1./Data_S(3,Data_S(7,:)==n_p);
% z = Data_S(8,Data_S(7,:)==n_p);
% % plot effects xi, for a given n. 
% % Select data
% i_xi  = unique(z); 
% markers_s = ["square","o","^","x"];
% colors_s  = ["#7E2F8E","#EDB120","#A2142F","#77AC30"];
% figure(1)
% hold on 
% for i = 1:length(i_xi)
%     s=scatter(x(z==i_xi(i)),y(z==i_xi(i)),20);
%     s.MarkerEdgeColor = 'k';
%     s.Marker          = markers_s{i};
%     s.MarkerFaceColor = colors_s{i};
%     s.LineWidth       = 0.2; 
% end
% legend({'$ \xi^{S} = 1.0$','$ \xi^{S} = 10.0 $','$ \xi^{S} = 100.0 $','$\xi^{S} = 1000.0$'},'Interpreter','latex');
% grid on
% box on
% set(gca, 'XScale', 'log')
% xlim([10^(-8),5])
% ylim([-0.1,1.1])
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
% print('../FigureS3_A','-dpng')
% clf;
% 
% x = [];
% y = [];
% z = []; 
% 
% xi  = 10.0; 
% x = Data_S(1,Data_S(8,:)==xi);
% y = Data_S(7,Data_S(8,:)==xi)./Data_S(3,Data_S(8,:)==xi);
% z = Data_S(7,Data_S(8,:)==xi);
% i_xi  = unique(z); 
% 
% figure(2)
% hold on 
% for i = 1:length(i_xi)
%     s=scatter(x(z==i_xi(i)),y(z==i_xi(i)),20);
%     s.MarkerEdgeColor = 'k';
%     s.Marker          = markers_s{i};
%     s.MarkerFaceColor = colors_s{i};
%     s.LineWidth       = 0.2; 
% end
% legend({(['$ n  =$',num2str(i_xi(1))]),(['$ n  =$',num2str(i_xi(2))]),(['$ n  =$',num2str(i_xi(3))])},'Interpreter','latex');
% grid on
% box on
% xlim([10^(-8),5])
% ylim([-0.1,7.5])
% set(gca, 'XScale', 'log')
% 
% %set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
% print('../FigureS3_B','-dpng')
% clf;
% 
% figure(3)
% y = 1./(((Data_S(3,Data_S(8,:)==xi))));
% 
% hold on 
% for i = 1:length(i_xi)
%     s=scatter(x(z==i_xi(i)),y(z==i_xi(i)),20);
%     s.MarkerEdgeColor = 'k';
%     s.Marker          = 'o';
%     s.MarkerFaceColor = colors_s{i};
%     s.LineWidth       = 0.2; 
% end
% legend({(['$ n  =$',num2str(i_xi(1))]),(['$ n  =$',num2str(i_xi(2))]),(['$ n  =$',num2str(i_xi(3))])},'Interpreter','latex');
% grid on
% box on
% xlim([10^(-8),5])
% ylim([-0.1,1.1])
% set(gca, 'XScale', 'log')
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
% print('../FigureS3_C','-dpng')
% 
% figure(5)
% x = Data_S(11,Data_S(8,:)==xi)./(365.25*60*60*1e6);
% hold on 
% for i = 1:length(i_xi)
%     s=scatter(x(z==i_xi(i)),(i_xi(i)./(y(z==i_xi(i))).*x(z==i_xi(i))),20);
%     s.MarkerEdgeColor = 'k';
%     s.Marker          = 'o';
%     s.MarkerFaceColor = colors_s{i};
%     s.LineWidth       = 0.2; 
% end
% legend({(['$ n  =$',num2str(i_xi(1))]),(['$ n  =$',num2str(i_xi(2))]),(['$ n  =$',num2str(i_xi(3))])},'Interpreter','latex');
% grid on
% box on
% %xlabel('$log_{10}(\Lambda)$[n.d.]',Interpreter='latex')
% %ylabel('$t_{det} [Myrs]$',Interpreter='latex')
% 
% %xlim([10^(-8),5])
% %ylim([0.3,10.1])
% %set(gca, 'XScale', 'log')
% %set(gca,'XTickLabel',[])
% %set(gca,'YTickLabel',[])
% print('../FigureS3_E','-dpng')
% 
