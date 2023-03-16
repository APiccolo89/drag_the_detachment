function [res,ID] = optimizing_fit(ID,nlm,f,D,t,tdet,type,test_Gif)
%=========================================================================$
%
%
%==========================================================================
ID.ID_A.fetch(1) = abs(f(1));
ID.ID_A.fetch(2) =(f(2));
[TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
if type == 1 
    [D,D_0D2D,t] = clean_data(TestData,D,t);
    res = goodnessOfFit(D,D_0D2D,'NMSE');
else
    if isempty(TestData.t_det)
        res = tdet;
    else
        res = tdet-TestData.t_det.*ID.n;
    end
end
    if test_Gif == 1

    figure(1)
    ax = gca; 
        hold on
        p(1) = plot(t.*ID.n,D,'LineStyle','-','Color','k','LineWidth',1.2); 
        if ID.ID_A.flag  == 1 
            p(3) = plot(t.*ID.n,D_0D2D,'LineStyle','-.','Color','r','LineWidth',1.0); 
        end
        p(2) = plot(t.*ID.n,D_0D2D,LineStyle=":",Color=[0.1 0.1 0.1],LineWidth=0.5);

    ax.YLim = [0.1,1.0];
    ax.YColor = [0.0 0.0 0.0];
    ax.XColor = [0.0 0.0 0.0]; 
    ax.LineWidth = 1.0;
    ax.Box = 'on'; 
    ax.XGrid = 'on';
    ax.YGrid = 'on'; 
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.XLabel.String = '$t^{\dagger}$';
    ax.YLabel.String = '$D^{\dagger}$'; 
    ax.YLabel.Interpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    
    if ID.ID_A.flag  == 1
    filename = (['Test_D',num2str(ID.ID_A.iteration)]); 
    folder   = '../Gif_Karls';
    if not(isfolder(folder))
        mkdir(folder);
    end
    pt=fullfile(folder,filename);
    set(gcf, 'Color', 'w')  
    print(pt,'-r600','-dpng')
    
    end

    
    figure(2)
    ax2 = gca;
    hold on
    scatter(res,f(1),'black','filled','d');
    ax2.XLabel.String= '$res$';
    ax2.YLabel.String= '$\omega_{\mathrm{depth}}$';
    ax2.YLabel.Interpreter = 'latex';
    ax2.XLabel.Interpreter = 'latex';
    ax2.XColor = [0.0 0.0 0.0]; 
    ax2.YColor = [0.0 0.0 0.0];    
    ax2.Box = 'on';
    ax2.LineWidth = 1.0; 
    ax2.XGrid = 'on';
    ax2.YGrid  = ' on'; 
    ax2.XMinorTick ='on';
    ax2.YMinorTick ='on'; 

    if ID.ID_A.flag  == 1
    filename = (['Test_F2',num2str(ID.ID_A.iteration)]); 
    folder   = '../Gif_Karls';
    if not(isfolder(folder))
        mkdir(folder);
    end
    pt=fullfile(folder,filename);
    set(gcf, 'Color', 'w')  
    print(pt,'-r600','-dpng')
    
    end

    
    figure(3)
    ax3 = gca;
    hold on
    scatter(res,f(2),'black','filled','o');
    ax3.XLabel.String= '$res$';
    ax3.YLabel.String= '$\omega_{s}$';
    ax3.YLabel.Interpreter = 'latex';
    ax3.XLabel.Interpreter = 'latex';
    ax3.XColor = [0.0 0.0 0.0]; 
    ax3.YColor = [0.0 0.0 0.0]; 

    ax3.Box = 'on';
    ax3.LineWidth = 1.0; 
    ax3.XGrid = 'on';
    ax3.YGrid  = ' on'; 
    ax3.XMinorTick ='on';
    ax3.YMinorTick ='on'; 
    

    if ID.ID_A.flag  == 1
    filename = (['Test_F1',num2str(ID.ID_A.iteration)]); 
    folder   = '../Gif_Karls';
    if not(isfolder(folder))
        mkdir(folder);
    end
    pt=fullfile(folder,filename);
    set(gcf, 'Color', 'w')  
    print(pt,'-r600','-dpng')
    end
    end



end
