
function [D,D_0D2D,t] = clean_data(TestData,D,t)

% Clean the input data: When I compile the data from the numerical
% simulation, i neglect the effects of the algorithm behind the detection
% of detachment and it appears that time to time within the selection area
% a piece of slab reappeared screwing up a bit the detection at high tc. 

ind_first_nan = find(isnan(D),1);
D = D(1:ind_first_nan);
t = t(1:ind_first_nan);
D(ind_first_nan)=0.1;

time_0D = real(TestData.time);
D_0D    = real(TestData.D_norm);
D_0D2D = interp1(time_0D,D_0D,t,'linear','extrap');
D_0D2D(D_0D2D<0.1)=0.0;

end