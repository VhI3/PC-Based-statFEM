function [ccy, pk1] = pk2cauchy( pk2, f )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Transform second Piola-Kirchhoff stress to cauchy stress
%  formula: ccy = fsf'/j
%  Input:
%    pk2 - second Piola-Kirchhoff stress
%    f -  deformation tensor
%  Output:
%    ccy - Cauchy stress
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%
% Project: PC-Based-statFEM 
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)


pk2 = [pk2(1), pk2(3);
    pk2(3), pk2(2)];
% cauchy = f*pk2*f'/det(f);
ccy = f*pk2*f';
ccy = voigt(ccy);
% 1st Piola-Kirchhoff
pk1 = pk2*f';
pk1 = voigt(pk1);
end


function V = voigt( M )
if size(M,1) == 2
    V = [M(1,1) M(2,2) M(1,2)]';
elseif size(M,1) == 3
    V = [M(1,1) M(2,2) M(3,3) M(1,2) M(2,3) M(1,3)]';
end

end
