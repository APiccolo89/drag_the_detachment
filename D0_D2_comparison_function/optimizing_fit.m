function [res] = optimizing_fit(ID,nlm,f,D,t)
%=========================================================================$
%
%
%==========================================================================
ID.ID_A.fetch(1) = f(1);
ID.ID_A.fetch(2) = f(2);
[TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
[D,D_0D2D] = clean_data(TestData,D,t);
res = goodnessOfFit(D,D_0D2D,'NMSE');
end
