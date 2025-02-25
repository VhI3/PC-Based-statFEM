function [displacement] = FEM_1DBar_Tipload_Solver(preProc_Variables, E_vector)
% FEM_1DBar_Tipload_Solver
% --------------------------------------------------------------
% Solves the linear elastic boundary value problem for a 1D tension bar
% under a tip load using the Finite Element Method (FEM).
%
% The function computes the displacement field by:
%   - Assembling the global stiffness matrix.
%   - Applying boundary conditions and external forces.
%   - Solving the system of linear equations for nodal displacements.
%
% Inputs:
%   preProc_Variables - Structure containing preprocessing parameters:
%       - GDof: Global degrees of freedom.
%       - numberElements: Number of finite elements.
%       - elementNodes: Element connectivity matrix.
%       - xi, w: Gauss integration points and weights.
%       - f_bar: Magnitude of the tip load.
%       - activeDofs: Active degrees of freedom (excluding fixed nodes).
%       - prescribedDofs: Degrees of freedom with prescribed displacements.
%       - A: Cross-sectional area.
%       - nodeCoordinates: Nodal coordinates along the bar.
%   E_vector - Vector of Young's modulus values for each node.
%
% Outputs:
%   displacement - Computed displacement vector for the entire bar.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Extract variables from input structure
GDof = preProc_Variables.GDof;
numberElements = preProc_Variables.numberElements;
elementNodes = preProc_Variables.elementNodes;
xi = preProc_Variables.xi;                      % Gauss integration points
w = preProc_Variables.w;                        % Gauss weights
f_bar = preProc_Variables.f_bar;                % Applied tip load
activeDofs = preProc_Variables.activeDofs;      % Active degrees of freedom
prescribedDofs = preProc_Variables.prescribedDofs; % Prescribed (fixed) degrees of freedom
A = preProc_Variables.A;                        % Cross-sectional area
nodeCoordinates = preProc_Variables.nodeCoordinates; % Nodal coordinates

%% Define shape functions and their derivatives
% Two-node linear Lagrangian shape functions (1D elements)
Nshape = @(xi) [-(xi - 1) / 2, (xi + 1) / 2];    % Shape functions N1, N2
Nshape_derivate = @(xi) [-1/2, 1/2];             % Derivatives dN1/dxi, dN2/dxi

%% Initialize global matrices
stiffness = zeros(GDof, GDof);                   % Global stiffness matrix
force = zeros(GDof, 1);                          % Global force vector

%% Assembly of the global stiffness matrix
for e = 1:numberElements
    % Retrieve nodal indices and corresponding degrees of freedom
    indice = elementNodes(e, :);                 % Nodes of the current element
    elementDof = indice;                         % Degrees of freedom of the element
    
    % Young's modulus for the element (linear interpolation)
    E_element = E_vector(indice);
    
    % Initialize element stiffness matrix
    ke = zeros(length(elementDof));
    
    % Gaussian quadrature for numerical integration
    for i = 1:length(xi)
        % Shape functions and interpolation of Young's modulus
        N1 = Nshape(xi(i));                      % Shape functions at Gauss point
        E_xi = N1(1)*E_element(1) + N1(2)*E_element(2); % Interpolated Young's modulus
        
        % Compute shape function derivatives and Jacobian transformation
        naturalDerivatives = Nshape_derivate(xi(i));
        [JacobianMatrix, ~, XYDerivatives] = Jacobian1D(nodeCoordinates(indice), naturalDerivatives);
        
        % Strain-displacement matrix (B matrix)
        B = XYDerivatives;                      % Derivative of shape functions in global coordinates
        
        % Element stiffness matrix contribution
        ke = ke + w(i) * (B' * E_xi * A * B) * det(JacobianMatrix);
    end
    
    % Assemble element stiffness into the global stiffness matrix
    stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + ke;
end

%% Application of external forces
force(end) = f_bar;                              % Apply tip load at the last node

%% Apply boundary conditions and solve
% Partition the global stiffness matrix and force vector
Krr = stiffness(activeDofs, activeDofs);         % Reduced stiffness matrix (active DOFs)
Kru = stiffness(activeDofs, prescribedDofs);    % Coupling between active and prescribed DOFs
Rr = force(activeDofs);                          % Reduced force vector

% Initialize displacement vector
displacement = zeros(GDof, 1);
Uu = displacement(prescribedDofs);               % Prescribed displacement (usually zero)

% Solve for active displacements
Ur = full(Krr \ (Rr - Kru * Uu));                % Solve the reduced system
displacement(activeDofs) = Ur;                   % Update displacement vector

end
