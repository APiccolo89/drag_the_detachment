classdef line_plot_post_process
    %LINE_PLOT_POST_PROCESS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        figure_number %number figure
        logx          %log scale
        logy          %log scale
        logcolor      %log color
        colormap_f      % colormap
        xlabel        %label x
        ylabel        %label y
        clabel        %label colorbar
        size_picture   %size figure
        save_path     %path save
        legend_option %legend if needed
        x             % x vector
        y             % y vector
        c             % z vector
        name_figure   % name figure
        clim          %limit of axis
        ctick         %ticks of the colorbar
        xlim          % lim x
        ylim          %lim y
    end

    methods
        function multiline_plot(obj,ntests,color_lists)
            figure(obj.figure_number)
            clf;

            set(gcf, 'Units','centimeters', 'Position', [0, 0, obj.size_picture(1),obj.size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [obj.size_picture(1), obj.size_picture(2)])
            ax = gca;
            for i = 1:ntests
                hold on
                aa = obj.x(:,i);
                bb = obj.y(:,i);
                cc = obj.c(:,i);
                aa = aa(~isnan(aa));
                bb = bb(~isnan(aa));
                cc = cc(~isnan(aa));
                if ~isempty(aa)
                    p2(i)=plot(aa,bb,'Color',obj.colormap_f(color_lists(i),:),'LineWidth',0.9);
                end
            end
            colormap(obj.colormap_f);
            ax.TickLabelInterpreter = 'latex';
            ax.LineWidth = 1.2;
            ax.Box = 'on';
            ax.FontSize = 14;
            % X Axis
            ax.XColor = [0,0,0];
            ax.XLabel.Interpreter = 'latex';
            % Y Axis
            ax.YColor =[0,0,0];
            ax.XColor = [0,0,0];
            ax.YLabel.Interpreter = 'latex';
            ax.XScale = obj.logx;
            ax.YScale = obj.logy;
            if ~isempty(obj.xlim)
                ax.XLim = [obj.xlim];
            end
            if ~isempty(obj.ylim)
                ax.YLim = [obj.ylim];
            end

            % Labels,Scales
            ax.XLabel.String = obj.xlabel;
            ax.YLabel.String = obj.ylabel;
            cb    = colorbar(gca,'southoutside');
            cb.Label.Interpreter = 'latex';
            cb.Label.String = obj.clabel;
            cb.Limits = [obj.clim];
            Ticks =  [obj.ctick];
            for i=1:numel(obj.ctick)
               Tickslabel{i} = num2str(obj.ctick(i));
            end
            cb.Ticks = Ticks;
            cb.TickLabelInterpreter = 'latex';
            cb.TickLabels = Tickslabel;
            cb.FontSize = 12;
            cb.Color    = [0,0,0];
            ax.Layer = 'top';
            caxis(obj.clim);
            filename = fullfile(obj.save_path);
            if ~isdir(filename)
                mkdir(filename);
            end
            filename = fullfile(filename,strcat(obj.name_figure,'.png'));
            exportgraphics(ax,filename,'Resolution',600,'BackgroundColor','white')
        end
    end
end

