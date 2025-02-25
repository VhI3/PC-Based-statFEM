function [H_s] = Hermite_noisePC(n_y, p_order)
% HERMIT_NOISEPC - Generates 1D Hermite polynomial basis functions for noise modeling.
%
% This function computes the Hermite polynomials using a recurrence relation and 
% generates polynomial chaos basis functions for a given number of random variables (RVs).
%
% Inputs:
%   n_y     - Number of random variables (dimensions).
%   p_order - Order of the Hermite polynomial chaos expansion.
%
% Outputs:
%   H_s     - Cell array of size (n_y x p_order) containing symbolic Hermite polynomials 
%             for each random variable.
%
% Notes:
%   - Hermite polynomials are orthogonal with respect to the standard Gaussian measure.
%   - The recurrence relation used is:
%       He_0(x) = 1  
%       He_1(x) = x  
%       He_{n+1}(x) = x * He_n(x) - (n - 1) * He_{n-1}(x)
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Step 1: Initialize and Compute 1D Hermite Polynomials using Recurrence
syms zetta;                     % Symbolic variable for polynomial generation
He_s = cell(p_order, 1);        % Cell to store Hermite polynomials
He_s{1} = sym(1);               % He_0(x) = 1
He_s{2} = zetta;                % He_1(x) = zetta

% Generate higher-order Hermite polynomials using recurrence:
for j = 2:p_order
    % Recurrence relation: He_{j+1}(x) = x * He_j(x) - (j-1) * He_{j-1}(x)
    He_s{j+1} = expand(zetta * He_s{j} - (j - 1) * He_s{j-1});
end

%% Step 2: Generate Hermite Polynomials for Each Random Variable
M = 1;                          % Number of random variables per dimension (fixed as 1)
H_s = cell(n_y, p_order);       % Initialize output cell array for polynomials

% For each random variable dimension:
for j = 1:n_y
    % Define a unique symbolic variable for the j-th random variable
    x{j} = sym(sprintf('zetta_%d', j));  
    
    % Assign the first-order Hermite polynomial to the j-th variable:
    % Note: Only the second polynomial (He_s{2}) is used here. 
    % To extend to higher orders, replace He_s{2} with appropriate He_s{p_order}.
    H_s{j} = subs(He_s{2}, zetta, x{j});
end

end
