%% plateWithHole_FEM_processLE.m
% --------------------------------------------------------------
% Perform Linear Elastic (LE) FEM Analysis for a Plate with a Hole.
%
% Description:
% This function computes the displacement field under linear elastic 
% conditions for two cases:
%   1. Homogeneous Young's modulus across the plate.
%   2. Inhomogeneous Young's modulus near the hole (to simulate punch effects).
%
% Inputs:
%   BVP : Structure containing preprocessed mesh, material, and boundary conditions.
%
% Outputs:
%   BVP : Updated structure with displacement results and related fields.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

function BVP = plateWithHole_FEM_processLE(BVP)

%% 1. Assign Variables from BVP Structure
elementNodes = BVP.preProc.msh.elementNodes;      % Element connectivity matrix
E_vector = BVP.preProc.material.E_vector;         % Homogeneous Young's modulus
E_punch_vector = BVP.preProc.material.E_punch_vector; % Inhomogeneous Young's modulus near the hole

%% 2. FEM Calculation (Homogeneous Material)
% Solve FEM problem assuming homogeneous Young’s modulus
[u, force] = FEM_plateWithHole_Solver(BVP.preProc, E_vector);

% Extract displacement components
ux = u(1:2:end, :);  % X-direction displacement
uy = u(2:2:end, :);  % Y-direction displacement

% Interpolate displacement values for surface plotting
ux_Surf = makeSurf(elementNodes, ux);
uy_Surf = makeSurf(elementNodes, uy);

%% 3. FEM Calculation (Inhomogeneous Material Near Hole)
% Solve FEM problem with varying Young’s modulus around the hole
[u_nonHom, force_nonHom] = FEM_plateWithHole_Solver(BVP.preProc, E_punch_vector);

% Extract displacement components for inhomogeneous case
ux_nonHom = u_nonHom(1:2:end, :);  % X-direction displacement
uy_nonHom = u_nonHom(2:2:end, :);  % Y-direction displacement

%% Assign back to BVP
BVP.proc.LE.u                           = u;
BVP.proc.LE.ux                          = ux;
BVP.proc.LE.uy                          = uy;
BVP.proc.LE.u_nonHom                    = u_nonHom;
BVP.proc.LE.force_nonHom                = force_nonHom;
BVP.proc.LE.ux_nonHom                   = ux_nonHom;
BVP.proc.LE.uy_nonHom                   = uy_nonHom;
BVP.proc.LE.ux_Surf                     = ux_Surf;
BVP.proc.LE.uy_Surf                     = uy_Surf;
BVP.proc.LE.force                       = force;

disp('2. Linear elastic FEM analysis for the Plate with a Hole is completed.');
end