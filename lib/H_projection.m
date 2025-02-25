function H = H_projection(senCoordinates, nodeCoordinates, elementNodes)
% This function computes the projection matrix H that maps the global displacement field
% to the displacement field at the sensor locations using linear interpolation.

% Inputs:
% senCoordinates: Coordinates of sensors
% nodeCoordinates: Coordinates of all nodes in the mesh
% elementNodes: Connectivity matrix defining which nodes belong to each element

% Outputs:
% H: Projection matrix

% Initialize the projection matrix H with zeros
% Project: PC-Based-statFEM 
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)


%
H = zeros(length(senCoordinates), size(nodeCoordinates, 1));

% Precompute element coordinates from the node coordinates and element connectivity
for i = 1:size(elementNodes, 1)
    elementCoordinates(i, :) = nodeCoordinates(elementNodes(i, :));
end

gridCell = [];

% Iterate over each sensor coordinate to determine its position relative to the mesh elements
for i = 1:length(senCoordinates)
    for j = 1:size(elementNodes, 1)
        % Check if the sensor is within the current element or exactly at the element's nodes
        if (senCoordinates(i) > elementCoordinates(j, 1) && senCoordinates(i) < elementCoordinates(j, 2)) ...
                || (senCoordinates(i) == elementCoordinates(j, 1) || (senCoordinates(i) == elementCoordinates(j, 2)))
            % If the sensor is within the element or at its boundary, store the relevant information
            gridCell = [gridCell; elementCoordinates(j, :) senCoordinates(i) elementNodes(j, :)];
        end
    end
end

% Compute the interpolation weights for each sensor based on its local coordinate within the element
for i = 1:size(gridCell, 1)
    % Transform the sensor's global coordinate into its local coordinate within the element
    xi = (2 * gridCell(i, 3) - gridCell(i, 1) - gridCell(i, 2)) / (gridCell(i, 2) - gridCell(i, 1));
    % Compute the shape functions for linear interpolation
    N1 = Lagrange1D(xi, 1, [-1 1]);
    % Assign the interpolation weights to the projection matrix
    H(i, gridCell(i, 4)) = N1(1);
    H(i, gridCell(i, 5)) = N1(2);
end

end
