%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clf



eta0 = 1e21; 
eta0D = 10*eta0;
L0    = 500e3;
D0    = 80e3; 
n     = 3.5; 
tau0  = 100e6;
drho   = 2*tau0/(9.81*L0);
B_n   = tau0^(1-n)/eta0;
B_d   = 1/eta0D; 
etaUM_vec =10.^(17:0.1:20); 
name_1 = num2str(log10(eta0));

D = D0;

for i = 1:length(etaUM_vec)
    etaUM = etaUM_vec(i); 
    Testdata = Compute_dddt_Drag(D,B_n,B_d,n,L0,D0,drho,etaUM,tau0,name_1);
    name_2 = append('T',name_1,num2str(i));
    Tests.(name_2) = Testdata;
    Testdata = [];
    eta_UM = [];  
end
i = 1 ; 
fn = fieldnames(Tests);
cc = jet(length(etaUM_vec));
for k = 1:numel(fn) 

TD = Tests.(fn{k});
T = TD(1,:);
D = TD(2,:);
hold on 
plot(T,D,'Color',cc(i,:))

set(gca, 'YScale', 'log')
grid on 
xlim([0,10])
ylim([10^(-2),10^(0)])

xlabel('t/tc [n.d.]')

i = i+1; 
end
print('Global_Test','-dpng')
clf; 
close; 



bla = 0.0