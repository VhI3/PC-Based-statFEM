%% plateWithHole_H_matrix.m
% --------------------------------------------------------------
% Statistical FEM: Projection Matrix (H) Computation
%
% Description:
% This function computes the projection matrix H that maps the global displacement
% field to sensor locations. It is essential for comparing simulation results
% with sensor measurements. The matrix H allows for extracting displacements at
% specific sensor nodes from the global solution vector.
%
% Inputs:
%   BVP : Boundary Value Problem structure with preprocessed FEM and sensor data.
%
% Outputs:
%   BVP : Updated structure containing:
%           - Full projection matrix (H_full)
%           - Projection matrices for active DOFs (Hx, Hy)
%           - Node indices corresponding to sensors (senNode)
%           - Reduced projection matrix (H_full_active)
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% -------------------------------------------------------------

function BVP = plateWithHole_H_matrix(BVP)
% This function computes the projection matrix H that maps the global displacement field to the sensor locations
% and applies this projection to the mean displacement field to obtain the projected mean displacement at sensor locations.

%% Assign from BVP
senCoor          = BVP.preProc.statfem.sensorCoordinates;
nSen             = BVP.preProc.statfem.nuSensors;
nodeCoordinates  = BVP.preProc.msh.nodeCoordinates;
GDOFs            = BVP.preProc.msh.GDOFs;
DOFs             = BVP.preProc.msh.DOFs;
activeDOFsX      = BVP.preProc.BC.activeDOFsX;
activeDOFsY      = BVP.preProc.BC.activeDOFsY;
activeDOFs       = BVP.preProc.BC.activeDOFs;
mu_u_pc_active = BVP.ssfem.mu_u_pc_active; % Mean displacement at active degrees of freedom

%% Calcualte
H_full = zeros(nSen*DOFs,GDOFs);
senNode = zeros(nSen,1);
cc = 1;
for i = 1:nSen
    rows = find(abs(nodeCoordinates(:,1)- senCoor(i,1)) < eps);
    nodes = abs(nodeCoordinates(rows,2) -  senCoor(i,2)) < eps;
    senNode(i) = rows(nodes);
    H_full(cc,DOFs*senNode(i)-1) = 1;
    H_full(cc+1,DOFs*senNode(i)) = 1;
    cc = cc +2;
end
Hx              = H_full(1:2:end,activeDOFsX);
Hy              = H_full(2:2:end,activeDOFsY);
Hx_full              = H_full(1:2:end,1:2:end);
Hy_full              = H_full(2:2:end,2:2:end);
P{1}          = Hx_full;
P{2}          = Hy_full;
H_full_active = H_full(:,activeDOFs);


%% Assign back to BVP
BVP.obs.senNode = senNode;
BVP.obs.H_full = H_full;
BVP.obs.Hx_full = Hx_full;
BVP.obs.Hy_full = Hy_full;
BVP.obs.Hx = Hx;
BVP.obs.Hy = Hy;
BVP.obs.P = P;
BVP.obs.H_full_active = H_full_active;
end