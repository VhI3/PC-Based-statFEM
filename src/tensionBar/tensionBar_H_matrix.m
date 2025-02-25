function BVP = tensionBar_H_matrix(BVP)
% tensionBar_H_matrix - Computes the projection matrix H to map the global displacement field
% to the sensor locations for the tension bar problem.
%
% The projection matrix H is essential in applications where measurements (e.g., sensor data)
% need to be compared or fused with the numerical model (FEM displacement field).
%
% This function:
%   - Calculates the full projection matrix H_full between all nodes and sensor locations.
%   - Extracts the active DOF projection matrix H (excluding the fixed boundary node).
%
% Inputs:
%   BVP - Boundary Value Problem structure containing:
%         * nodeCoordinates: Node positions in the FEM model.
%         * elementNodes: Element connectivity (node pairs per element).
%         * senCoordinates: Sensor positions along the bar.
%
% Outputs (updated in BVP structure):
%   BVP.statfem.H        - Projection matrix for active DOFs (size: nSensors x nActiveDOFs).
%   BVP.statfem.H_full   - Full projection matrix including fixed DOFs (size: nSensors x nNodes).
%
% Note:
%   * The first node is assumed to be fixed (zero displacement), so it's excluded from H.
%   * The function `H_projection` (called internally) performs the actual construction of H.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% -------------------------------------------------------------------------

%% Extracting Inputs from BVP
nodeCoordinates = BVP.preProc.nodeCoordinates;   % FEM nodal coordinates [nNodes x 1]
elementNodes    = BVP.preProc.elementNodes;      % Element-node connectivity [nElements x 2]
senCoordinates  = BVP.statfem.senCoordinates;    % Sensor coordinates [nSensors x 1]

%% Projection Matrix Calculation
% H_full: Maps global displacement vector (including fixed nodes) to sensor displacements.
% Each row of H_full corresponds to a sensor location; columns correspond to FEM nodes.
H_full = H_projection(senCoordinates, nodeCoordinates, elementNodes);

% Extract projection matrix for active degrees of freedom.
% Since the first node is fixed, we remove its contribution (first column).
H = H_full(:, 2:end); 

%% Update BVP Structure with Results
BVP.statfem.H      = H;       % Projection matrix for active DOFs
BVP.statfem.H_full = H_full;  % Full projection matrix (including fixed node)

end