function BVP = tensionBar_FEM_process_LE(BVP, calcase)
% tensionBar_FEM_process_LE
% --------------------------------------------------------------
% Processes the linear elastic boundary value problem (BVP) for a
% one-dimensional tension bar under a tip load using the Finite Element Method (FEM).
%
% This function:
%   - Extracts geometry, material properties, and boundary conditions from the BVP structure.
%   - Solves the displacement field using a linear FEM approach.
%   - Computes the analytical solution for comparison.
%   - Stores the computed solutions back into the BVP structure.
%
% Inputs:
%   BVP - Boundary Value Problem structure containing preprocessing data.
%
% Outputs:
%   BVP - Updated structure containing:
%         - FEM displacement solution
%         - Analytical solution (homogeneous displacement field)
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Extract necessary variables from BVP structure
preProc_Variables = BVP.preProc;                  % Extract preprocessing data
mu_E = preProc_Variables.mu_E;                    % Mean Young's modulus [GPa]
nodeCoordinates = preProc_Variables.nodeCoordinates; % Node coordinates [mm]
f_bar = preProc_Variables.f_bar;                  % Applied tip load [kN]
A = preProc_Variables.A;                          % Cross-sectional area [mm²]
activeDofs = preProc_Variables.activeDofs;        % Active degrees of freedom

%% Perform FEM calculation
% Solve the linear elastic FEM problem for the tension bar under a tip load.
switch calcase
    case 'linear'
        E_vector = mu_E * ones(preProc_Variables.GDof, 1); % Young's modulus vector for all elements
    case 'nonlinear'
        E_func = @(x) mu_E*((sin(x/10)/1.5) + 1);
        E_vector  = E_func(nodeCoordinates);
    otherwise
        error('Unknown calculation case type selected.');
end

% FEM_1DBar_Tipload_Solver should return the global displacement vector.
U_FEM = FEM_1DBar_Tipload_Solver(preProc_Variables, E_vector');

%% Compute analytical solution
% Homogeneous solution:
% Formula: u(x) = (f_bar * x) / (mu_E * A)
U_analytic = (f_bar * nodeCoordinates) / (mu_E * A);

% Extract displacements at active degrees of freedom
U_analytic_active = U_analytic(activeDofs);

%% Assign results back to the BVP structure
BVP.fem.proc.displacement = U_FEM;                % FEM displacement solution
BVP.analyticalSolution.U_analytic = U_analytic;   % Full analytical solution
BVP.analyticalSolution.U_analytic_active = U_analytic_active; % Active DOFs analytical solution

disp('2. Finished processing the linear elastic boundary value problem for a tension bar using FEM.');
end
