

function [II]= plot_2Dmaps_profiles(S,pt_save,nlm)
close all;
clf;
FIT = S.FIT;
TB  = S.TB;
font_axes = 16;
font_legend = 14;
font_text   = 5;
size_shit = [12,13.5];
LineWidth = 1.2;
marker_size = 10;
fnames = fieldnames(TB);
z = -1000:1050/512:50;
z = z(z<20);
ind_z = find(z>-200);
flength = length(fnames);

for ktest = 1:flength
    if TB.(fnames{ktest}).P_Var.failed == 0
        t = TB.(fnames{ktest}).t.*TB.(fnames{ktest}).ID.n;
        D_m = TB.(fnames{ktest}).D_matrix;
        D_m(D_m==-Inf)=nan;
        for i=1:length(t)
            a = D_m(i,:);
            D_mv(i,:) = movmean(a,10);
            b = D_mv(i,:);
            ind=find(b(ind_z(1):ind_z(end))==nanmin(b(ind_z(1):ind_z(end))),1);
            x(i) = z(ind_z(ind(1)));
            DD(i)=b(ind_z(ind(1)));
        end
        if strcmp(fnames{ktest},'T22')
            bla=0;
        end
%         figure(1)
%         clf;
%         ax = gca;
%         pcolor(t,z,D_m(1:length(t),:)'./80);shading flat; colorbar;
%         colorbar;
%         caxis([0.1,1.0])
%         colormap(crameri('nuuk',27))
%         hold on
%         ylim([-300,-100])
%         ax.Box = 'on';
%         ax.LineWidth = 1.2;
%         ax.XColor = [0 0 0 ];
%         ax.YColor = [0 0 0 ];
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.Layer = 'top';
%         figure_name = fnames{ktest};
%         pt_fig = fullfile(pt_save,figure_name);
%         print(pt_fig,'-dpng','-r600');
%         bla =0;
% 
%         figure(2)
%         clf;
%         ax = gca;
%         hold on
% 
%         plot(D_mv./80,z,'k')
%         plot(DD./80,x,'r','LineWidth',1.2)
%         yline(TB.(fnames{ktest}).Dp2-100,Color='b',LineStyle=':',LineWidth=1.2)
%         yline(-100-TB.(fnames{ktest}).f(1).*TB.(fnames{ktest}).P_Var.L0./1e3,Color='r',LineStyle='--',LineWidth=1.2);
%         hold on
%         yline(TB.(fnames{ktest}).Dp2-100,Color='b',LineStyle=':',LineWidth=1.2)
%         yline(-100-TB.(fnames{ktest}).f(1).*TB.(fnames{ktest}).P_Var.L0./1e3,Color='r',LineStyle='--',LineWidth=1.2)
%         xlim([0.01,0.9])
%         ylim([-300,-100])
%         ax.Box = 'on';
%         ax.LineWidth = 1.2;
%         ax.XColor = [0 0 0 ];
%         ax.YColor = [0 0 0 ];
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         figure_name = strcat(figure_name,'profile');
%         pt_fig = fullfile(pt_save,figure_name);
%         print(pt_fig,'-dpng','-r600');
% 
%         figure(3)
%         hold on
        buf1 = -100-TB.(fnames{ktest}).f(1).*TB.(fnames{ktest}).P_Var.L0./1e3;
         buf2 = TB.(fnames{ktest}).Dp2-100;
%         I = (x-buf1)./(mean(x(end-1))-buf1);
%         I = movmean(I,4);
%         plot(t,log10(abs(I)))
        II(ktest)=(buf1-mean(x(end-1)))./(mean(x(end-1)));

        x = [];
        DD= [];
        D_mv = [];
        I = [];

    end
end


end
