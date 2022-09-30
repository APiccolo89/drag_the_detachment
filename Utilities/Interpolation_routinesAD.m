function [Interpolation] = Interpolation_routinesAD(Testdata,Testdata_a,Benchmark)
% 
% Input 
%==========================================================================
% Testdata = Dimensional results 
% Testdata_a = Dimensional results
%==========================================================================
% Output
% Interpolation vectors
%==========================================================================
    D_D = Testdata.D_norm;
    t_D = Testdata.time;
    D_A = Testdata_a.D_norm;
    t_A = Testdata_a.time; 
    D_Di = interp1(t_D,D_D,t_A);
    D_Ai = interp1(t_A,D_A,t_D);
    error_DA = nanmean(abs(D_D-D_Ai));
    error_AD = nanmean(abs(D_Di-D_A));
    if (Benchmark==1)
        disp([':::::::::::::::::::::::::::::::::::::::::::::::::::'])
        disp(['Errors interpolation are: ']);
        disp(['Adimensional interpolated to dimensional:',num2str(error_DA,'%3e')])
        disp(['Dimensional interpolated to adimensional:',num2str(error_AD, '%3e')])
        disp(['Error between detected detachment times:',num2str(abs(Testdata_a.t_det-Testdata.t_det),'%3e')])
        disp(['=================================================='])
    end
    Interpolation.t_A = t_A; 
    Interpolation.t_D = t_D;
    Interpolation.D_Di = D_Di; 
    Interpolation.D_Ai = D_Ai; 
    Interpolation.error_DA = [error_DA,nanmax(abs(D_D-D_Ai)),nanmin(abs(D_D-D_Ai))]; 
    Interpolation.error_AD = [error_AD,nanmax(abs(D_Di-D_A)),nanmin(abs(D_Di-D_A))]; 
end
