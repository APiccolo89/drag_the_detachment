classdef maps_plot
    %MAPS_PLOT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        figure_number %number figure

        logx          %log scale

        logy          %log scale

        logcolor      %log color

        colormap_f    % colormap

        xlabel        %label x

        ylabel        %label y

        clabel        %label colorbar

        size_picture   %size figure

        save_path     %path save

        legend_option %legend if needed

        x             % x vector

        y             % y vector

        c             % z vector

        contourf_option      % use contourf options

        clim          %limit of the colorbar

        name_figure   % name figure
    end

    methods
        function make_maps(obj)

            figure(obj.figure_number)

            clf;

            set(gcf, 'Units','centimeters', 'Position', ...
                [0, 0, obj.size_picture(1),obj.size_picture(2)], ...
                'PaperUnits', 'centimeters', 'PaperSize', ...
                [obj.size_picture(1), obj.size_picture(2)])

            ax = gca;

            if obj.logcolor == 1

                obj.c = log10(obj.c);
            end

            if obj.contourf_option{1}==1

                 if ~isempty(obj.contourf_option{2})
                     levels = obj.clim(1):(obj.clim(2)-obj.clim(1))/obj.contourf_option{2}:obj.clim(2);
                 end
                p1 = contourf(obj.x,obj.y,obj.c,obj.contourf_option{2});



            else

                p1 = pcolor(obj.x,obj.y,obj.c);

                shading interp;
            end

            caxis(obj.clim);

            ax.XScale = obj.logx;

            ax.YScale = obj.logy;

            ax.XLabel.Interpreter = 'latex';

            ax.YLabel.Interpreter = 'latex';

            ax.XAxis.TickLabelInterpreter = 'latex';

            ax.YAxis.TickLabelInterpreter = 'latex';

            ax.XLabel.String = obj.xlabel;

            ax.YLabel.String = obj.ylabel;

            ax.LineWidth = 1.2;

            ax.FontSize = 12;

            ax.Box = 'on';

            ax.Layer = 'top';

            ax.XColor = [0,0,0];

            ax.YColor = [0,0,0];

            ax.TickDir = 'both';

            path2colormap = strcat('Utilities\ScientificColourMaps8\',obj.colormap_f,'\',obj.colormap_f,'.mat');

            load(path2colormap);

            colormap(eval(obj.colormap_f));

            cb = colorbar('southoutside');

            cb.Label.Interpreter = 'latex';

            cb.Label.String = obj.clabel;

            cb.TickLabelInterpreter = 'latex';

            cb.FontSize = 10;

            filename = fullfile(obj.save_path);

            if ~isfolder(filename)

                mkdir(filename);
            end

            filename = fullfile(filename,strcat(obj.name_figure,'.png'));

            exportgraphics(ax,filename,'Resolution',600,'BackgroundColor' ...
                ,'white')

        end

    end
end

