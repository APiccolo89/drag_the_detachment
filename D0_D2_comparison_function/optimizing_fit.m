function [res,ID] = optimizing_fit(ID,nlm,f,D,t,tdet,type,test_Gif,number_fetches,matver)
%=========================================================================$
%
%
%==========================================================================
for i = 1:number_fetches
    ID.ID_A.fetch(i) = abs(f(i));

end
res = [];
while isempty(res)
    try
        [TestData] = Run_Simulation_DragA(ID.ID_A,nlm);
        if type == 1
            [D,D_0D2D,t] = clean_data(TestData,D,t);
            z_norm_real = (D-mean(D))./std(D);
            z_norm_test = (D_0D2D-mean(D_0D2D))./std(D_0D2D);
            res = goodnessOfFit(D,D_0D2D,'NRMSE');
            if matver<2020

                res = abs(1-res);
            else

                res = abs(res);
            end


            
        else
            if isempty(TestData.t_det)
                res = tdet;
            else
                res = tdet-TestData.t_det.*ID.n;
            end
            [D,D_0D2D,t] = clean_data(TestData,D,t);
            res2 = goodnessOfFit(D,D_0D2D,'NRMSE');

        end
    catch
        bla=0
    end




end
