%%
function [Table,L,td_rel_E]=Create_Table_Latex(TB,RES,fetch,pt_save)
% fieldname 
field_names = fieldnames(TB);
for i = 1:length(field_names)
    T = TB.(field_names{i});
    % Find Lambda
    Lambda{i} = sprintf('%0.2f',log10(T.ID.Lambda));
    L(i) =log10(T.ID.Lambda); 
    % Find xi um
    if isnan(T.P_Var.xiUM)
        xiUM{i} = sprintf('%0.2f',0.0);
        if RES(i)>0
         R{i} = '$HR$' ;
        else
          R{i} ='$LR$';
        end
    else
        xiUM{i} =sprintf('%0.2f',T.P_Var.xiUM);
        R{i} = '$HR$'; 
    end
    % Find xi S
    xiS{i} = sprintf('%0.2f',10.0);
    %find eta0DS
    eta0DS{i} = sprintf('%0.2f',log10(T.P_Var.eta0DS));
    %find eta0DM
    eta0DM{i} = sprintf('%0.2f',log10(T.P_Var.eta0DM));
    %find tau0
    tau0{i} = sprintf('%0.2f',(T.P_Var.s0)./1e6);
    %find l0
    L0{i} = sprintf('%0.2f',T.P_Var.L0./1e3);
    %find tc
    tc{i} = sprintf('%0.2f',T.P_Var.tc);
    %find td
    td_buf = T.Detachment_0D;
    if isempty(td_buf)
        td_2D{i} = '$[]$';
        td_0D{i} = '$[]$';
        td_rel_err{i} = '$[]$';
        Res{i} = '$[]$';
        if fetch==2 
                    fit_par1{i} = '$[]$';
                                fit_par2{i} = '$[]$';

        else
        fit_par{i} = '$[]$';
        end
        td_rel_E(i) =nan;
    else
        td_2D{i} = sprintf('%0.2f',T.Detachment_2D);
        td_0D{i} = sprintf('%0.2f',T.Detachment_0D);
        td_rel_err{i} =sprintf('%0.2f', abs((T.Detachment_2D-T.Detachment_0D)./(T.Detachment_2D))*100);
        td_rel_E(i) = abs((T.Detachment_2D-T.Detachment_0D)./(T.Detachment_2D))*100; 
        if fetch ==2 
        fit_par1{i} = sprintf('%0.2f',T.f(1));
        fit_par2{i} = sprintf('%0.2f',T.f(2));
        else
        fit_par{i} = sprintf('%0.2f',T.f);
        end
        Res{i} = sprintf('%0.2f',T.res*100);
    end
end
if fetch ==2 
        Table = table(td_rel_err',fit_par1',fit_par2',Res');

else
    Table = table(Lambda',xiUM',xiS',eta0DM',eta0DS',tau0',L0',tc',td_2D',td_0D',td_rel_err',fit_par',Res',R');
end
Table.Properties.RowNames = field_names; 
% Use the 

table2latex(Table, pt_save)


end

function table2latex(T, filename)
% ----------------------------------------------------------------------- %
% Function table2latex(T, filename) converts a given MATLAB(R) table into %
% a plain .tex file with LaTeX formatting.                                %
%                                                                         %
%   Input parameters:                                                     %
%       - T:        MATLAB(R) table. The table should contain only the    %
%                   following data types: numeric, boolean, char or string.
%                   Avoid including structs or cells.                     %
%       - filename: (Optional) Output path, including the name of the file.
%                   If not specified, the table will be stored in a       %
%                   './table.tex' file.                                   %  
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       LastName = {'Sanchez';'Johnson';'Li';'Diaz';'Brown'};             %
%       Age = [38;43;38;40;49];                                           %
%       Smoker = logical([1;0;1;0;1]);                                    %
%       Height = [71;69;64;67;64];                                        %
%       Weight = [176;163;131;133;119];                                   %
%       T = table(Age,Smoker,Height,Weight);                              %
%       T.Properties.RowNames = LastName;                                 %
%       table2latex(T);                                                   %                                       
% ----------------------------------------------------------------------- %
%   Version: 1.1                                                          %
%   Author:  Victor Martinez Cagigal                                      %
%   Date:    09/10/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
    
    % Error detection and default parameters
    if nargin < 2
        filename = 'table.tex';
        fprintf('Output path is not defined. The table will be written in %s.\n', filename); 
    elseif ~ischar(filename)
        error('The output file name must be a string.');
    else
        if ~strcmp(filename(end-3:end), '.tex')
            filename = [filename '.tex'];
        end
    end
    if nargin < 1, error('Not enough parameters.'); end
    if ~istable(T), error('Input must be a table.'); end
    
    % Parameters
    n_col = size(T,2);
    col_spec = [];
    for c = 1:n_col, col_spec = [col_spec 'l']; end
    col_names = strjoin(T.Properties.VariableNames, ' & ');
    row_names = T.Properties.RowNames;
    if ~isempty(row_names)
        col_spec = ['l' col_spec]; 
        col_names = ['& ' col_names];
    end
    
    % Writing header
    fileID = fopen(filename, 'w');
    fprintf(fileID, '\\begin{tabular}{%s}\n', col_spec);
    fprintf(fileID, '%s \\\\ \n', col_names);
    fprintf(fileID, '\\hline \n');
    
    % Writing the data
    try
        for row = 1:size(T,1)
            temp{1,n_col} = [];
            for col = 1:n_col
                value = T{row,col};
                if isstruct(value), error('Table must not contain structs.'); end
                while iscell(value), value = value{1,1}; end
                if isinf(value), value = '$\infty$'; end
                temp{1,col} = num2str(value);
            end
            if ~isempty(row_names)
                temp = [row_names{row}, temp];
            end
            fprintf(fileID, '%s \\\\ \n', strjoin(temp, ' & '));
            clear temp;
        end
    catch
        error('Unknown error. Make sure that table only contains chars, strings or numeric values.');
    end
    
    % Closing the file
    fprintf(fileID, '\\hline \n');
    fprintf(fileID, '\\end{tabular}');
    fclose(fileID);
end



