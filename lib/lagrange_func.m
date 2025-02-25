function [ N ] = lagrange_func( xi, xik, p )
% This function produces the value of Lagrange Shape functions.
% Input : xi ==============> The point on which the derivative shape
%                            function has to be calculated.
%       : xik==============> The range of xi. Normaly between [-1 1].
%       : p  ==============> Degree of polynomials.
% Output:   Matrix of shape functions, N(i,j) is the value of 
%           the i-th shape function at postition xi(j)
% Author: Vahab B.Narouie
% Last Change : 18. February 2019
%

%%
% check the compatible
if(size(xik,2)~=(p+1))
    return;
end
%initial lize the shape functions matrix
N = zeros(p+1,size(xi,2));
%loop over all the xi values
for j=1:size(xi,2)
    %loop over the shape functions
    for i=1:p+1
        N(i,j)=1;
        %loop k
        for k=1:p+1
            %and skip k=i
            if(k ~= i)
                N(i,j)=N(i,j)*(xik(k)-xi(j))/(xik(k)-xik(i));
                %fprintf("i=%d \t, k=%d \t xik(k)=%f \t xi=%f \t xik(i)=%f \t \n",i,k,xik(k),xi(j),xik(i));
            end
        end
    end
end
%
end

