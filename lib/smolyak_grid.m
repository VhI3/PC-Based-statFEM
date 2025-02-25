function [xd, wd] = smolyak_grid(d, stages, oned_rule_func)
% SMOLYAK_GRID - Constructs a Smolyak sparse grid for multidimensional quadrature.
%
% Syntax:
%   [xd, wd] = smolyak_grid(d, stages, oned_rule_func)
%
% Inputs:
%   d             - Number of dimensions.
%   stages        - Number of refinement levels (accuracy level).
%   oned_rule_func - Function handle or cell array of function handles for 1D quadrature rule.
%
% Outputs:
%   xd - Matrix of size (d x n) containing Smolyak grid points.
%   wd - Column vector of length n containing the corresponding weights.
%
% Example:
%   [xd, wd] = smolyak_grid(2, 3, @gauss_hermite_rule);
%
% Description:
%   This function constructs a Smolyak sparse grid using the given one-dimensional
%   quadrature rule function(s). The Smolyak grid is built from tensor products of 
%   1D quadrature rules and weighted combinations to achieve a more efficient 
%   quadrature for high-dimensional problems.
%
% Project: PC-Based-statFEM
% Forker: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% Originally created: %   Elmar Zander
% Institute of Scientific Computing, TU Braunschweig.
% --------------------------------------------------------------

%% Ensure oned_rule_func is a cell array with the correct size
if ~iscell(oned_rule_func)
    oned_rule_func = repmat({oned_rule_func}, d, 1);
elseif length(oned_rule_func) == 1
    oned_rule_func = repmat(oned_rule_func, d, 1);
elseif length(oned_rule_func) ~= d
    error('Dimension d does not match the size of the cell array of quadrature functions.');
end

%% Generate 1D quadrature formulas for each dimension and level
x1 = cell(d, stages); % Store quadrature nodes
w1 = cell(d, stages); % Store quadrature weights
n1d = zeros(d, stages); % Number of points in each quadrature rule

for k = 1:d
    for j = 1:stages
        [x1{k, j}, w1{k, j}] = funcall(oned_rule_func{k}, j);
        n1d(k, j) = length(x1{k, j});
    end
end

%% Construct Smolyak grid nodes and weights
xd = [];
wd = [];

% Generate multi-index sets for Smolyak quadrature
I = multiindex(d, stages - 1, 'combine', false);
I = cell2mat(I(max(0, stages - d) + 1:end)') + 1;

for i = 1:size(I, 1)
    alpha = I(i, :); % Multi-index set

    % Select the corresponding 1D quadrature nodes and weights
    tmp_x1 = cell(d, 1);
    tmp_w1 = cell(d, 1);
    
    for j = 1:d
        alpha_j = alpha(j);
        tmp_x1{j} = x1{j, alpha_j};
        tmp_w1{j} = w1{j, alpha_j};
    end
    
    % Compute the tensor product of selected quadrature rules
    [tmp_xd, tmp_wd] = tensor_mesh(tmp_x1, tmp_w1);

    % Compute the Smolyak combination factor
    alpha_sum = sum(alpha);
    factor = (-1)^(d + stages - 1 - alpha_sum) * nchoosek(d - 1, alpha_sum - stages);

    % Append new points and weights
    xd = [xd, tmp_xd]; %#ok<AGROW>
    wd = [wd; factor * tmp_wd]; %#ok<AGROW>
end

%% Remove duplicate points and sum the corresponding weights
[xdt, ~, j] = unique(xd', 'rows'); % Identify unique grid points
xd = xdt';
wd = accumarray(j, wd); % Sum weights for duplicate points

end
