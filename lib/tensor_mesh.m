function [xd, wd] = tensor_mesh(x1, w1)
% TENSOR_MESH - Generates a D-dimensional tensor product grid from 1D meshes and weights.
%
% Syntax:
%   [xd, wd] = tensor_mesh(x1, w1)
%   xd = tensor_mesh(x1)
%
% Inputs:
%   x1 - Cell array containing 1D quadrature nodes for each dimension.
%   w1 - (Optional) Cell array containing 1D quadrature weights.
%
% Outputs:
%   xd - Matrix of size (d x nd) containing the tensor product nodes.
%   wd - Column vector of length nd containing the corresponding weights (if provided).
%
% Example:
%   [x1, w1] = gauss_hermite_rule(10);
%   [x2, w2] = gauss_hermite_rule(12);
%   [xd, wd] = tensor_mesh({x1, x2}, {w1, w2});
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

% Determine if weights should be computed
compute_weights = (nargin == 2);

% Input validation
d = numel(x1);  % Number of dimensions
if compute_weights && d ~= numel(w1)
    error('Dimension mismatch: x1 and w1 must have the same number of elements.');
end

% Compute the number of points in each dimension
points_per_dim = cellfun(@numel, x1);  

% Compute total number of tensor product points
total_points = prod(points_per_dim);

% Preallocate the tensor grid and weights (if applicable)
xd = zeros(d, total_points);
if compute_weights
    wd = ones(total_points, 1);
end

% Generate tensor product grid using efficient index expansion
grid_idx = cell(1, d);
[grid_idx{:}] = ndgrid(x1{:});  % Generate multi-dimensional grid

% Flatten the grid into a matrix format
for k = 1:d
    xd(k, :) = grid_idx{k}(:).';
end

% Compute weights if required
if compute_weights
    weight_idx = cell(1, d);
    [weight_idx{:}] = ndgrid(w1{:});  % Compute weight grid
    wd = prod(cell2mat(cellfun(@(w) w(:), weight_idx, 'UniformOutput', false)), 2);
end

end
