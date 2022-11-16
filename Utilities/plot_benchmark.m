function  plot_benchmark(Intp,ID,D,Dim)
    pt0 ='../Utilities/benchmark_results';
    if(Dim == 1)
    t_ = Intp.t_D;
    D_ = Intp.D_Ai;
    pt=fullfile(pt0,'Benchmark_A_D.png');
    pt2=fullfile(pt0,'Benchmark_A_D_diff.png');

    labels={'Dimensional','Adimensional','Analytical'};
    name = strcat("$log_{10}(\Lambda) = $",num2str(log10(ID.Lambda),3));
    

    else
    t_ = Intp.t_A;
    D_ = Intp.D_Di;  
    pt=fullfile(pt0,'Benchmark_D_A.png');
    pt2=fullfile(pt0,'Benchmark_D_A_diff.png');

    labels={'Adimensional','Dimensional','Analytical'};
    name = strcat("Benchmark, $\Lambda = $",num2str(log10(ID.Lambda),3), ' [n.d.], Dimensional to Adimensional');

    end
    if not(isfolder(pt0))
        mkdir(pt0)

    end
    

    figure(1)
    % Collect the field names 
    plot(ID.n.*t_,D,"Color",'k','LineStyle','--','LineWidth',1.5)
    hold on 
    plot(ID.n.*t_,D_,"Color",'red','LineStyle','-','LineWidth',1.0)
    % Analytical solution
    t=(0.0:1e-5:1/ID.n);
    D_anal = (1-ID.n.*t).^(1/ID.n);
    plot(ID.n.*t,D_anal,"Color",'blue','LineStyle','-','LineWidth',0.8)
    legend(labels)
    grid on
    %grid minor
    xlabel('$\frac{t}{t_d} [n.d.]$',Interpreter='latex')
    ylabel('$\frac{D}{D_0} [n.d.]$',Interpreter='latex')
    xlim([0,10])
    ylim([0.1,1.0])
    title(name,Interpreter="latex")
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'Color','none')
    print(pt,'-dpng')
    clf; 
    close;

    figure(1)
    % Collect the field names 
    plot(t_,(abs(D_-D)),"Color",'red','LineStyle','-','LineWidth',1.0)
    grid on
    grid minor
    if Dim == 1 
        ylabel('$ abs|\left(\frac{D}{D_0}\right)^{A}-\left(\frac{D}{D_0}\right)^{D}| [n.d.]$',Interpreter='latex')

    else
        ylabel('$ abs|\left(\frac{D}{D_0}\right)^{D}-\left(\frac{D}{D_0}\right)^{A}| [n.d.]$',Interpreter='latex')

    end
    xlabel('$\frac{t}{t_c} [n.d.]$',Interpreter='latex')
    set(gca, 'YScale', 'log')
    title(name,Interpreter="latex")
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'Color','none')
    print(pt2,'-dpng')
    clf; 
    close;

end