classdef scatter_plot_post_process
    properties
        figure_number %number figure
        logx          %log scale
        logy          %log scale
        logcolor      %log color
        colormap_f      % colormap
        colormap_discrete % number of element
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
        xlim 
        ylim 
        loop_unique_c
        clim 
        size 
        markers 
        conditions 
    end

    methods
        function make_scatter_plot(obj)
        if isempty(obj.size)
            s = 60;
        else
            s = obj.size; 
        end

        figure(obj.figure_number)
        clf; 
        set(gcf, 'Units','centimeters', 'Position', [0, 0, obj.size_picture(1),obj.size_picture(2)], 'PaperUnits', 'centimeters', 'PaperSize', [obj.size_picture(1), obj.size_picture(2)])

        ax = gca;
        if length(obj.c) == length(obj.x)
            if obj.logcolor
                obj.c = log10(obj.c);
            end
        else 
            disp('The figure will not have any colorbar');
        end
        if ~isempty(obj.loop_unique_c)
            val = unique(obj.c);
          for i = 1:length(val)
            color = obj.colormap_f(i,:); 
            hold on
            p1 = scatter(obj.x(obj.c==val(i)),obj.y(obj.c==val(i)),s,color,'filled','MarkerEdgeColor','k');
          end
          hold off
        elseif (~isempty(obj.markers))
            size_marker = numel(obj.markers);
            condition_vector = obj.conditions{1};
            condition_condition = obj.conditions{2}; 
            allowed = condition_vector<=condition_condition;

            
             p1 = scatter(obj.x(allowed==1),obj.y(allowed==1),s,obj.c(allowed==1),'filled','MarkerEdgeColor','k','Marker','hexagram')
                hold on 
                p2 = scatter(obj.x(allowed==0),obj.y(allowed==0),s,obj.c(allowed==0),'filled','MarkerEdgeColor','k','Marker','d')
        
        
        else
          p1 = scatter(obj.x,obj.y,s,obj.c,'filled','MarkerEdgeColor','k');

        end
        ax.XScale = obj.logx;
        ax.YScale = obj.logy;
        if ~isempty(obj.colormap_f) && isempty(obj.loop_unique_c)
            % load colormap 
            path2colormap = strcat('Utilities\ScientificColourMaps8\',obj.colormap_f,'\',obj.colormap_f,'.mat'); 
            load(path2colormap);
            values = eval(obj.colormap_f);
            if ~isempty(obj.colormap_discrete)
                P = size(values,1);
                values = interp1(1:size(values,1), values, linspace(1,P,obj.colormap_discrete), 'linear');
            end
            colormap(values)
            cb = colorbar();
            cb.Label.Interpreter = 'latex';
            cb.Label.String = obj.clabel;
            cb.TickLabelInterpreter = 'latex';
            cb.FontSize = 10;
            if ~isempty(obj.clim)
                ax.CLim = obj.clim; 
            end
        end
        ax.Box = 'on';
        ax.LineWidth = 1.2;
        ax.XColor = [0,0,0];
        ax.YColor = [0,0,0];
        ax.XLabel.Interpreter = 'latex';
        ax.YLabel.Interpreter = 'latex';
        ax.XLabel.String = obj.xlabel;
        ax.YLabel.String = obj.ylabel;
        ax.XAxis.TickLabelInterpreter = 'latex';
        ax.YAxis.TickLabelInterpreter = 'latex';
        ax.FontSize = 12;
        if ~isempty(obj.xlim)
            ax.XLim = obj.xlim; 
        end
        if ~isempty(obj.ylim)
            ax.YLim = obj.ylim;
        end
        filename = fullfile(obj.save_path);
        if ~isdir(filename)
            mkdir(filename);
        end
        filename = fullfile(filename,strcat(obj.name_figure,'.png'));
        exportgraphics(ax,filename,'Resolution',600,'BackgroundColor','white')         
        end
    end
end

