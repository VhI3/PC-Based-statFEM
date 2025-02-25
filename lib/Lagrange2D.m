function [ N ] = Lagrange2D( xi1,xi2,p1,p2,xik1,xik2 )
%% This function produces the value of Lagrange Shape functions.
% Input : xi1 ==============> The point on which the shape function has to 
%                             be calculated.in X (or X1).
%       : xi2 ==============> The point on which the shape function has to
%                             be calculated.in Y (or X2).
%       : xik1==============> The range of xi1. Normaly between [-1 1].
%       : xik2==============> The range of xi2. Normaly between [-1 1].
%       : p1  ==============> Degree of polynomials. in X (or X1).
%       : p2  ==============> Degree of polynomials. in Y (or X2).
% Output:   Vector of shape functions postition xi1(i), xi2(j)
% FOR 2-D Problems.
%
% Author: Vahab B.Narouie
% Last Change : 18. February 2019
%

%%
%       p1=1,p2=1
%
%   3------------4                 4------------3
%   |            |      order      |            |
%   |            |     ===>        |            |
%   |            |                 |            |
%   1------------2                 1------------2
%
%       p1=1,p2=2
%      
%   5------------6                4------------3
%   |            |     order      |            |
%   3            4     ===>       6            5 
%   |            |                |            |
%   1------------2                1------------2
%
%       p1=2,p2=2
%      
%   7-----8------9                4-----7------3 
%   |            |    order       |            |
%   4     5      6    ===>        8     9      6
%   |            |                |            |
%   1-----2------3                1-----5------2
%
%       p1=3,p2=2
%      
%   9--10---11---12               4---9----8---3
%   |            |      order     |            | 
%   5   6    7   8     ===>      10   11   12  7
%   |            |                |            |
%   1---2----3---4                1---5----6---2
%
%
%
%       p1=3,p2=3
%      
%   13-----14-----15------16              4-----10-------9------3
%   |                     |               |                     |
%   9      10     11      12     order   11     16       15     8
%   |                     |      ===>     |                     |
%   5      6      7       8              12     13       14     7
%   |                     |               |                     |
%   1------2------3-------4               1------5-------6------2
%%
N = zeros(1,(p1+1)*(p2+1));
N1 = lagrange_func(xi1,xik1,p1);
N2 = lagrange_func(xi2,xik2,p2);
c = 1;
for j=1:p2+1
    for i=1:p1+1
        N(1,c) = N1(i,1)*N2(j,1);
        c = c +1;
    end
end
%
end


%% TODO
%  This function is checked Till p1 = 1, p2 = 1  with [-1 1] and [-1 1].
%  1. Check it for p1 = 2,3 ,  p2 = 2,3 with [-1 1] and [-1 1].
%  2. Check it for p1 = 1, p2 = 1  with [0 1] and [0 1].
%  3. Check it for p1 = 1, p2 = 1  with [-1 1] and [0 1].
%  4. Check it for p1 = 2, p2 = 1  with [-1 1] and [0 1].


