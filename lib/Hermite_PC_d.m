function [H_s_meanFree,c_ij] = Hermite_PC_d(M,p_order)
% Project: PC-Based-statFEM 
% Forker: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% Originally created by Felipe Uribe (Oct/2014)
% Changes and enhancements: Added comments, cleaned code structure,
%                           updated formatting, and improved readability.


%% Calculate 1D Hermite polynomials: Recurrence relation
% symbolic
syms chi;
He_s    = cell(p_order,1);
He_s{1} = sym(1);
He_s{2} = chi;
for j = 2:p_order
    He_s{j+1} = expand(chi*He_s{j} - (j-1)*He_s{j-1});
end
% Define the number of RVs
x   = cell(1,M);
H_s = cell(p_order,M);   % Hermite polinomial for each dimension
for j = 1:M
    x{j} = sym(sprintf('chi_%d',j));
    for i = 1:p_order+1
        H_s{i,j} = subs(He_s{i},chi,x{j});
    end
end
H_s_meanFree = H_s(2:end, :);

H_ss = reshape(H_s_meanFree,size(H_s_meanFree,1)*size(H_s_meanFree,2),1);
H_s0 = [H_ss H_ss];
% %%
% %% M-dimensional PC computation
Psi_s  = cell(size(H_s0,1)^2,1);   % symbolic version
cc = 1;
for i = 1:size(H_s0,1)
    for j = 1:size(H_s0,1)
        Psi_s{cc} = H_s0{i,1}*H_s0{j,2};
        cc = cc +1;
    end
end
x = cell(1,M);
for j = 1:M
    x{j} = sym(sprintf('chi_%d',j));
end
%
dw   = @(x) exp(-0.5*x^2)/sqrt(2*pi);  %% Gaussian measure
P = p_order*M;
c_ij = zeros(P,P);
end
