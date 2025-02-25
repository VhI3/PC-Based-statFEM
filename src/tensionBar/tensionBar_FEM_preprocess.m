function BVP = tensionBar_FEM_preprocess(BVP, i_le, calcase)
%% tensionBar_FEM_preprocess
% --------------------------------------------------------------
% Pre-processing function for the 1D tension bar under a tip load problem.
%
% This function sets up the geometry, mesh, material properties, and
% random field models for uncertainty quantification via Karhunen-Loève (KL)
% expansions. The uncertain Young's modulus is modeled as a weakly-stationary
% Gaussian random field with a specified correlation length and covariance kernel.
%
% Inputs:
%   BVP  - Boundary Value Problem structure (empty or partially filled)
%   i_le - Index specifying the correlation length for the Young's modulus
%
% Outputs:
%   BVP  - Updated structure with preprocessing information, including:
%          - Geometry and mesh properties
%          - Material uncertainty and random field parameters
%          - KL expansion parameters
%          - Observation noise properties
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Geometry Properties
DIM = 1;           % Problem dimension (1D)
L = 100;           % Length of the bar [mm]
A = 20;            % Cross-sectional area [mm²]
f_bar = 800;       % Distributed force along the bar [kN]

%% Mesh Properties
numberElements = 99;                              % Number of elements
numberNodes = numberElements + 1;                 % Total number of nodes
DOFs = 1;                                         % Degrees of freedom per node
GDof = DOFs * numberNodes;                        % Global degrees of freedom
n_u = numberNodes - 1;                            % Number of active DOFs

nodeCoordinates = linspace(0, L, numberNodes)';   % Node coordinates
elementNodes = [(1:numberElements)', (2:numberElements+1)'];  % Element connectivity
xi_Range = [-1 1];                                % Parametric element domain
[prescribedDofs, activeDofs] = deal(1, setdiff(1:GDof, 1));  % Boundary conditions

% Gauss quadrature
numberGausPoints = 1;
[xi, w] = Gauss_int(numberGausPoints);

%% Uncertain Material Properties (Young's Modulus)
mu_E = 200;             % Mean Young's modulus [GPa]
sig_E = 0.2 * mu_E;     % Standard deviation [GPa]
mu_kappa = log(mu_E^2 / sqrt(mu_E^2 + sig_E^2));    % Log-normal mean
sig_kappa = sqrt(log(1 + (sig_E^2 / mu_E^2)));      % Log-normal std dev

%% KL Expansion of Young's Modulus
corrLength_E_vector = [10, 25, 50, 100];              % Correlation length options [mm]
corrLength_E = corrLength_E_vector(i_le);             % Selected correlation length
dom_bound = {[0, L]};                                  % Domain boundaries
p_order = 2;                                          % Polynomial order for KL

%% Kernel (Covariance) Functions
% The covariance function models the spatial correlation of the random field.
% Choose from the following kernels by setting `kernel_type`:

% kernel_type options:
% 'gaussian'  - Gaussian (squared exponential) kernel (infinitely differentiable)
% 'matern_5_2' - Matérn kernel with ν = 5/2 (smooth but less than Gaussian)
% 'matern_3_2' - Matérn kernel with ν = 3/2 (less smooth than 5/2)
% 'exponential' - Exponential kernel (less smooth, good for rough fields)

switch calcase
    case 'linear'
        kernel_type = 'gaussian'; % Select kernel type
    case 'nonlinear'
        kernel_type = 'matern_5_2'; % Select kernel type
    otherwise
        error('Unknown calculation case type selected.');
end

switch kernel_type
    case 'gaussian'
        % Gaussian Kernel (Squared Exponential):
        % cov(x1, x2) = exp(-0.5 * (|x1 - x2| / corrLength_E)^2)
        cov_func = @(x1, x2) exp(-0.5 * ((abs(x1 - x2)) / corrLength_E)^2);

    case 'matern_5_2'
        % Matérn Kernel (ν = 5/2):
        % cov(x1, x2) = (1 + sqrt(5)*r + 5r²/3) * exp(-sqrt(5)*r), where r = |x1 - x2|/corrLength_E
        cov_func = @(x1, x2) (1 + sqrt(5)*abs(x1 - x2)/corrLength_E + 5/3*(abs(x1 - x2)^2)/(corrLength_E^2)) ...
            .* exp(-sqrt(5)*abs(x1 - x2)/corrLength_E);

    case 'matern_3_2'
        % Matérn Kernel (ν = 3/2):
        % cov(x1, x2) = (1 + sqrt(3)*r) * exp(-sqrt(3)*r)
        cov_func = @(x1, x2) (1 + sqrt(3)*abs(x1 - x2)/corrLength_E) ...
            .* exp(-sqrt(3)*abs(x1 - x2)/corrLength_E);

    case 'exponential'
        % Exponential Kernel:
        % cov(x1, x2) = exp(-|x1 - x2| / corrLength_E)
        cov_func = @(x1, x2) exp(-abs(x1 - x2)/corrLength_E);

    otherwise
        error('Unknown kernel type selected.');
end

%% KL Expansion of Model-Reality Mismatch
P_d = 1;                  % Polynomial order for model-reality mismatch
nMC = 1e3;                % Number of Monte Carlo simulations

%% Observation Data Parameters
rho_preDefined = 1.5;     % Scaling factor for observations
switch calcase
    case 'linear'
        sig_e = 0.1;
    case 'nonlinear'
        sig_e = 1;
    otherwise
        error('Unknown calculation case type selected.');
end
% Standard deviation of observation noise
P_e = 1;                  % Polynomial order of noise

% Predefined eigenvalues for model-reality mismatch modes
tmp_eta_matrix_preDefined = [3; 3; 2.5; 2.3; 2.2; 0.7; 0.4; 0.3; 0.2; 0.1; 0.1; 0.05; 0.02; 0.01];

%% Assign Values to BVP Structure
% Geometry and loading
BVP.preProc.geometry = struct('DIM', DIM, 'L', L);
BVP.preProc.A = A;
BVP.preProc.f_bar = f_bar;

% Mesh information
BVP.preProc.numberElements = numberElements;
BVP.preProc.numberNodes = numberNodes;
BVP.preProc.DOFs = DOFs;
BVP.preProc.xi_Range = xi_Range;
BVP.preProc.n_u = n_u;
BVP.preProc.nodeCoordinates = nodeCoordinates;
BVP.preProc.GDof = GDof;
BVP.preProc.elementNodes = elementNodes;
BVP.preProc.numberGausPoints = numberGausPoints;
BVP.preProc.xi = xi;
BVP.preProc.w = w;
BVP.preProc.prescribedDofs = prescribedDofs;
BVP.preProc.activeDofs = activeDofs;

% Material uncertainty
BVP.preProc.mu_E = mu_E;
BVP.preProc.sig_E = sig_E;
BVP.preProc.mu_kappa = mu_kappa;
BVP.preProc.sig_kappa = sig_kappa;

% KL expansion of Young's modulus
BVP.preProc.corrLength_E = corrLength_E;
BVP.preProc.dom_bound = dom_bound;
BVP.preProc.p_order = p_order;
BVP.preProc.cov_func = cov_func;
BVP.preProc.kernel_type = kernel_type;
BVP.preProc.nMC = nMC;

% Model-reality mismatch
BVP.preProc.P_d = P_d;
BVP.preProc.rho_preDefined = rho_preDefined;
BVP.preProc.tmp_eta_matrix_preDefined = tmp_eta_matrix_preDefined;

% Observation noise
BVP.preProc.P_e = P_e;
BVP.preProc.sig_e = sig_e;

disp('1. Pre-processing for the Tension Bar problem is completed.');
end
