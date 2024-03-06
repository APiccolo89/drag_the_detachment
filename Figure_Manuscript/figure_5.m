function figure_5(Data_,ptsave)
% Create Data Set 
%=========================================================================%
% Test_2D -> D,tau
% Test fitting -> D,tau 
%=========================================================================%

for ktest = 1:flength
    if ~isempty(TB.(fnames{ktest}).D)
        t = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).ID.n;
        if nlm.Linear == 0
            Lambda = (TB.(fnames{ktest}).ID.Lambda)./(1+TB.(fnames{ktest}).ID.Df_UM.*TB.(fnames{ktest}).ID.tau_mc.^(3.5-1));
        else
            Lambda = (TB.(fnames{ktest}).ID.Lambda);
        end


        % Fill the vector with dimensionless arrays: [2D numerical test]
        %=====================================================================%
        F_S(1,1:length(t),ktest)= t;
        F_D(1,1:length(t),ktest)= t;
        F_Tnd(1,1:length(t),ktest)= t;
        F_dTdt_nd(1,1:length(t),ktest)= t;
        %=====================================================================%
        % Fill the vector with dimensional time
        t_dim = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).tc;
        t_dim = t_dim./(365.25*60*60*24*1e6);
        %=====================================================================%
        F_T(1,1:length(t),ktest)= t_dim;
        F_dTdt(1,1:length(t),ktest)= t_dim;
        %=====================================================================%
        % Fill the fitting time
        t0D = TB.(fnames{ktest}).t0D.*TB.(fnames{ktest}).ID.n;
        %======================================================================%
        F_Sft(1,1:length(t0D),ktest)= t0D;
        F_Dft(1,1:length(t),ktest)= t;
        %Fill the properties
        D = TB.(fnames{ktest}).D;
        D_0D = TB.(fnames{ktest}).D0D2;
        tau = TB.(fnames{ktest}).tau;
        tau0D = TB.(fnames{ktest}).tau0D;
        %=====================================================================%
        F_S(2,1:length(t),ktest) = tau(1:length(t));
        F_Sft(2,1:length(t0D),ktest) = tau0D(3,:);
        F_D(2,1:length(t),ktest) = D(1:length(t));
        F_Dft(2,1:length(t),ktest) = D_0D;
        %=====================================================================%
        %=====================================================================%
        %Fill the color
        F_S(3,1:length(t),ktest) =ones(length(t),1).*TB.(fnames{ktest}).res.*100;
        F_Sft(3,1:length(t0D),ktest) =ones(length(t0D),1).*TB.(fnames{ktest}).res.*100;
        F_D(3,1:length(t),ktest) =ones(length(t),1).*TB.(fnames{ktest}).res.*100;
        F_Dft(3,1:length(t),ktest) =ones(length(t),1).*TB.(fnames{ktest}).res.*100;

        %=====================================================================%

        topo = TB.(fnames{ktest}).topo;
        dt_dim = t_dim(2:1:end)-t_dim(1:1:end-1);
        dTopo = topo(2:1:length(t_dim))-topo(1:1:length(t_dim)-1);
        dU = dTopo./dt_dim;
        T_Mean = (t_dim(2:1:end)+t_dim(1:end-1)).*0.5;
        t_M = 0.5.*(t(1:1:end-1)+t(2:1:end));
        %=====================================================================%
        F_Tnd(1,1:length(t),ktest)= t;
        F_dTdt_nd(1,1:length(t_M),ktest)= t_M;
        F_Tnd(2,1:length(t),ktest)= topo(1:length(t))./(-min(topo(1:length(t))));
        F_dTdt_nd(2,1:length(t_M),ktest)= dU;
        %
        F_T(1,1:length(t_dim),ktest)= t_dim;
        F_dTdt(1,1:length(T_Mean),ktest)= T_Mean;
        F_T(2,1:length(t_dim),ktest)= topo(1:length(t_dim))./(-min(topo(1:length(t))));
        F_dTdt(2,1:length(T_Mean),ktest)= dU;
        %
        F_T(3,1:length(T_Mean),ktest)= log10(ones(length(T_Mean),1).*Lambda);
        F_dTdt(3,1:length(T_Mean),ktest)= log10(ones(length(T_Mean),1).*Lambda);
        F_Tnd(3,1:length(T_Mean),ktest)= log10(ones(length(T_Mean),1).*Lambda);
        F_dTdt_nd(3,1:length(T_Mean),ktest)= log10(ones(length(T_Mean),1).*Lambda);
    end
end
path2colormap = strcat('Utilities\ScientificColourMaps8\','lipari','\','lipari','.mat');
load(path2colormap);
cmap2 = colormap(lipari);

path2colormap2 = strcat('Utilities\ScientificColourMaps8\','glasgow','\','glasgow','.mat');
load(path2colormap2);
cmap = colormap(glasgow);

fun_0D = Manuscript_function_Container;  

c0D = squeeze(F_Dft(3,:,:));
c2D = squeeze(F_D(3,:,:));

color_lists_FIT0D = fun_0D.color_computation(flength,c0D,0,20,0);

figure(1)
subplot(2,2,1)
subplot(2,2,2)
subplot(2,2,3)
subplot(2,2,4)



end

