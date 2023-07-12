%=========================================================================%

clear all;
close all;
addpath Adimensional\
addpath Dimensional\
addpath \Users\'Andrea Piccolo'\Dropbox\freezeColors-master\
addpath Utilities\
addpath D0_D2_comparison_function\
clear all
% Common Data
D0=80e3;
n=3.5;
Df_S=10;
nlm = Problem_type.Linear;
Df_UM = nan;
DB_path = '../Data_Base/Test_DB2.hdf5';
pt_save = 'C:\Users\Andrea Piccolo\Dropbox\Bayreuth_Marcel\BasicCode\Viscous_DRAG\Manuscript\Linear_Comparison_fetch';
file_name_table = 'Linear_Experiments_Table'
pt_Table = fullfile(pt_save,file_name_table);
if not(isdir(pt_save))
    mkdir(pt_save);
end
l = h5info(DB_path,'/Viscous/HR');
%l2 =h5info(DB_path,'/Viscous/LR_');
TestsA  = {l.Groups.Name};
%TestsB  = {l2.Groups.Name};
HR = 1:length(TestsA);
%LR = -(1:length(TestsB));
Res = [HR];
Tests   = [TestsA];

[TB,FIT] = perform_optimization_DataBase(Tests,n,D0,Df_S,nlm,Df_UM,DB_path,pt_save,Res,1);

LN_ALT = struct("TB",TB,"FIT",FIT);
[Table,L,td_rel_E]=Create_Table_Latex(TB,Res,0,pt_Table);

save('Data_Base_Fit.mat',"LN_ALT",'-append');
