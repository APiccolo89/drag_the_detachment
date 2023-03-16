function [B] = correct_data(B, n, r,p,Water,d0)

B = B*2^(n-1);
B = B.*10^(-n.*6); 
B = B.*d0.^(-p); 
B = B.*Water^(r); 
end

