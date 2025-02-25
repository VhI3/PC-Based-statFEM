function [eigval, eigvec_f, eigvec, KLterms] = KL_fredholm_Galerkin_local_1D(preProc_Variables, cov_func)
% KL_fredholm_Galerkin_local_1D
% --------------------------------------------------------------
% Performs the Karhunen-Loève (KL) expansion by solving the Fredholm
% integral equation using the Galerkin finite element method in 1D.
%
% The KL expansion provides a reduced-order representation of a stochastic
% process, useful for modeling spatially varying random fields such as
% material properties (e.g., Young's modulus).
%
% Inputs:
%   preProc_Variables - Structure containing:
%       - GDof: Global degrees of freedom.
%       - numberElements: Number of elements in the mesh.
%       - elementNodes: Connectivity matrix (elements to nodes).
%       - nodeCoordinates: Coordinates of the mesh nodes.
%   cov_func - Function handle defining the covariance function of the
%              random process (e.g., Gaussian or exponential covariance).
%
% Outputs:
%   eigval     - Vector of eigenvalues obtained from the KL expansion.
%   eigvec_f   - Cell array of function handles for eigenfunctions (continuous form).
%   eigvec     - Matrix of eigenfunctions evaluated at the nodal coordinates.
%   KLterms    - Number of terms retained in the KL expansion based on the energy ratio criterion.
%
% Project: PC-Based-statFEM
% Forker: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% Originally created by Felipe Uribe (Oct/2014)
% Changes and enhancements: Added comments, cleaned code structure,
%                           updated formatting, and improved readability.
% --------------------------------------------------------------

%% Extract variables from preprocessing structure
GDof = preProc_Variables.GDof;
numberElements = preProc_Variables.numberElements;
elementNodes = preProc_Variables.elementNodes;
nodeCoordinates = preProc_Variables.nodeCoordinates;

% Energy ratio for KL truncation: Determines how much variance to capture
EnergyRatio = 0.99;

%% Shape functions and Gauss quadrature setup
% Linear two-node Lagrangian shape functions
Nshape = @(xi) [-(xi - 1) / 2, (xi + 1) / 2];     % Shape functions N1, N2
Nshape_derivate = @(xi) [-1/2, 1/2];              % Shape function derivatives

% Gauss-Legendre quadrature (for accurate integration)
nGp = 6;                                         % Quadrature order
[x_gl, w_gl] = Gauss_int(nGp);                   % Quadrature points and weights

%% Assembly of Mass Matrix (B)
% The mass matrix B arises from the inner product of shape functions.
B = zeros(GDof, GDof);                           % Initialize global mass matrix

for e = 1:numberElements
    indice = elementNodes(e, :);                 % Nodes of element e
    elementDof = indice;                         % Element DOFs
    Be = zeros(length(elementDof));              % Element mass matrix
    
    for i = 1:length(x_gl)
        N1 = Nshape(x_gl(i));                    % Shape functions at quadrature point
        naturalDerivatives = Nshape_derivate(x_gl(i)); % Derivatives
        [JacobianMatrix, ~, ~] = Jacobian1D(nodeCoordinates(indice), naturalDerivatives); % Jacobian
        
        % Element mass matrix contribution
        Be = Be + w_gl(i) * (N1' * N1) * det(JacobianMatrix);
    end
    
    % Assemble element mass matrix into the global mass matrix
    B(elementDof, elementDof) = B(elementDof, elementDof) + Be;
end

%% Assembly of Covariance Matrix (C)
% The covariance matrix C encodes the covariance structure of the random field.
C = zeros(GDof, GDof);

for e = 1:numberElements
    indice_e = elementNodes(e, :);               % Element e node indices
    for f = 1:numberElements
        indice_f = elementNodes(f, :);           % Element f node indices
        Cef = zeros(length(indice_e));           % Element covariance matrix
        
        for i = 1:length(x_gl)
            Ne = Nshape(x_gl(i));                % Shape functions (element e)
            [JacobianMatrix_e, ~, ~] = Jacobian1D(nodeCoordinates(indice_e), Nshape_derivate(x_gl(i)));
            xp1 = Ne * nodeCoordinates(indice_e); % Physical coordinate (element e)
            
            for j = 1:length(x_gl)
                Nf = Nshape(x_gl(j));            % Shape functions (element f)
                [JacobianMatrix_f, ~, ~] = Jacobian1D(nodeCoordinates(indice_f), Nshape_derivate(x_gl(j)));
                xp2 = Nf * nodeCoordinates(indice_f); % Physical coordinate (element f)
                
                % Compute covariance contribution
                Cef = Cef + w_gl(i) * w_gl(j) * cov_func(xp1, xp2) * ...
                    (Ne' * Nf) * det(JacobianMatrix_e) * det(JacobianMatrix_f);
            end
        end
        
        % Assemble into global covariance matrix
        C(indice_e, indice_f) = C(indice_e, indice_f) + Cef;
    end
end

%% Solve the generalized eigenvalue problem: C * A = B * A * D
% Here, eigenvalues (D) represent variances, and eigenvectors (A) represent spatial modes.
[~, Afull] = eig(C, B);                          % Solve full generalized eigenvalue problem
cumulated_eigs = cumsum(abs(diag(Afull))) / sum(abs(diag(Afull))); % Cumulative variance captured

% Determine the number of KL terms needed to reach the energy ratio
KLterms = find(cumulated_eigs >= EnergyRatio, 1);

% Solve for the leading KLterms eigenvalues and eigenvectors
[D, A] = eigs(sparse(C), sparse(B), KLterms);    % Efficient computation for large systems
[eigval, idx] = sort(diag(A), 'descend');        % Sort eigenvalues in descending order
eigvec = D(:, idx);                              % Corresponding eigenvectors

%% Generate continuous eigenfunction handles for interpolation
eigvec_f = cell(KLterms, 1);                     % Cell array for eigenfunction handles

for k = 1:KLterms
    eigvec_f{k} = @(xx) interp1(nodeCoordinates, eigvec(:, k), xx, 'linear', 'extrap');
end

end
