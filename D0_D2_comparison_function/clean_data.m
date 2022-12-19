
function [D,D_0D2D] = clean_data(TestData,D,t)

time_0D = TestData.time;
D_0D    = TestData.D_norm;
D_0D2D = interp1(time_0D,D_0D,t,'linear','extrap');
D_0D2D(isnan(D_0D2D) | D_0D2D<0.1)=0.0;
D(isnan(D) | D<0.1)=0.0;
end