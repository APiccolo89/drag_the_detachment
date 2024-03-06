%=========================================================================%

clear all;
close all;
addpath Adimensional\

addpath Dimensional\


addpath Utilities\

addpath D0_D2_comparison_function\


clear all
% Matlab version command: It appears that the fitting option where a mess
% before Matlab 2020, and basically the fitting of certain methods where
% computed as % fitting and not % of the complement error. This wonderful
% discrepancy allowed me to spend an entire weekend finding the reason why
% the code was not working, till, for a chance, I saw the wonderful warning
% in the matlab website - offcourse, was not evident, and offcourse was at
% the very end of the page. To avoid further problem, I introduce a version
% control such that the output is consistent with the current version. I
% wanted to be extremely lengthy just for documents this issue. 
matlab_version = version('-release'); 
matlab_version = str2num(matlab_version(1:4));


% Common Data
D0=80e3;

n=3.5;

Df_S=10;

nlm = Problem_type;

nlm.Linear=0;   % Switching the position of linear-non_linear activate the non linear upper mantle routine.

nlm.iteration = 1; % Iterate for the stress
 
nlm.cut_off   = 1; % Introduce cut-off for the upper mantle {must work on the slab}

Df_UM = nan;  % Assume xium = to nan

DB_path = '../Data_Base/Test_DB2.hdf5'; % Database path

pt_save = '..\Viscous_DRAG\Manuscript\New_Fit_Definitive\';

file_name_table = 'Table_Experiments';

pt_Table = fullfile(pt_save,file_name_table);

if not(isdir(pt_save))

    mkdir(pt_save);
end

pt_save = fullfile(pt_save,'Fetches_2');

if not(isdir(pt_save))

    mkdir(pt_save);
end


l = h5info(DB_path,'/Viscous/HR');
l2 = h5info(DB_path,'/Viscous/HR_DIS');

TestsA  = {l.Groups.Name};
TestsB = {l2.Groups.Name};

Dif = ones(length(TestsA),1);
Dis = zeros(length(TestsB),1);

Res = [];

Rheo = [Dif;Dis];

Tests   = [TestsA,TestsB];

[TB,FIT] = perform_optimization_DataBase(Tests,n,D0,Df_S,nlm,Df_UM,DB_path,pt_save,Res,2,matlab_version,Rheo);

NLN_fetches2 = struct('TB',TB,'FIT',FIT);

[Table,L,td_rel_E]=Create_Table_Latex(TB,Res,2,pt_Table);

save('Data_Base_Fit_DEF.mat','NLN_fetches2');
