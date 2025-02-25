function [ N ] = Lagrange1D( xi1,p1,xik1 )
% This function produces the value of Lagrange Shape functions.
% Input : xi1 ==============> The point on which the shape function has to 
%                             be calculated.
%       : xik1==============> The range of xi1. Normaly between [-1 1].
%       : p1  ==============> Degree of polynomials. in X (or X1).
% Output:   Vector of shape functions, is the value of 
%           the i-th shape function at postition xi1(j) 
% FOR 1-D Problems.
%
% Author: Vahab B.Narouie
% Last Change : 18. February 2019
%

%%
N = zeros(1,(p1+1));
N1 = lagrange_func(xi1,xik1,p1);
c = 1;
for i=1:p1+1
    N(1,c) = N1(i,1);
    c = c +1;
end

end

%% TODO
%  This function is checked Till p1 = 3 with [-1 1].
%  1. Check it for p1 = 4,5,6,7,8,9,10 with [-1 1].
%  2. Check it for p1 = 1,2,3 with [0 1].
