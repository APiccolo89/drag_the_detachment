function [res,ID] = optimizing_fit(ID,nlm,f,D,t,tdet,type,test_Gif,number_fetches)
%=========================================================================$
%
%
%==========================================================================
for i = 1:number_fetches
    ID.ID_A.fetch(i) = abs(f(i));

end

[TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
if type == 1 
    [D,D_0D2D,t] = clean_data(TestData,D,t);
    res = goodnessOfFit(D,D_0D2D,'NRMSE');
else
    if isempty(TestData.t_det)
        res = tdet;
    else
        res = tdet-TestData.t_det.*ID.n;
    end
end




end
