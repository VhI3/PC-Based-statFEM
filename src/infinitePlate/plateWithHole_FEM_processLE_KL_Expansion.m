%% plateWithHole_FEM_processLE_KL_Expansion.m
% --------------------------------------------------------------
% Karhunen-Loève Expansion & Polynomial Chaos for Young's Modulus Uncertainty
%
% Description:
% This function performs the Karhunen-Loève (KL) expansion to model the
% uncertainty in Young's modulus for a 2D plate with a hole.
% It uses the Fredholm-Galerkin method and a Polynomial Chaos (PC) expansion
% to represent the stochastic field with random variables.
%
% Steps:
% 1. Perform KL expansion to obtain eigenvalues and eigenvectors.
% 2. Calculate Polynomial Chaos basis functions and norms.
% 3. Generate Smolyak quadrature points and weights.
% 4. Generate stochastic field realizations of Young's modulus.
%
% Inputs:
%   BVP : Boundary Value Problem structure with preprocessing variables.
%
% Outputs:
%   BVP : Updated structure containing KL, PC, and stochastic field results.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

function BVP = plateWithHole_FEM_processLE_KL_Expansion(BVP)

%% 1. Extract Variables from BVP
preProc_Variables = BVP.preProc;
L = preProc_Variables.geometry.L;           % Plate length
H = preProc_Variables.geometry.H;           % Plate height
nuNodes = preProc_Variables.msh.nuNodes;    % Number of mesh nodes
mu_kappa = preProc_Variables.UQ.mu_kappa;   % Mean of the stochastic field
sig_kappa = preProc_Variables.UQ.sig_kappa; % Standard deviation of the stochastic field
p_order = 2;                                % Polynomial Chaos expansion order
Psi_com = [];                               % Initialize matrix for combined Psi values

%% 2. Karhunen-Loève Expansion Setup
% Set correlation lengths in x and y directions
Lx_kappa = L / 2;
Ly_kappa = H / 2;

% Define the covariance function using a Matérn 5/2 kernel
cov_func = @(x1, x2) ...
    (1 + sqrt(5)*abs(x1(1) - x2(1))/Lx_kappa + (5/3)*(abs(x1(1) - x2(1))^2)/(Lx_kappa^2)) * ...
    (1 + sqrt(5)*abs(x1(2) - x2(2))/Ly_kappa + (5/3)*(abs(x1(2) - x2(2))^2)/(Ly_kappa^2)) * ...
    exp(-sqrt(5)*(abs(x1(1) - x2(1))/Lx_kappa + abs(x1(2) - x2(2))/Ly_kappa));

% Perform KL expansion to obtain eigenvalues (eigval), eigenvectors (eigvec), and the number of terms (M_kappa)
[eigval, ~, eigvec, M_kappa] = KL_fredholm_Galerkin_local_2D_optimized(preProc_Variables, cov_func);

%% 3. Polynomial Chaos Expansion (PCE)
% Compute Hermite polynomial chaos basis functions and norms
[~, Psi_u, ~, PsiSqNorm, P_u] = Hermite_PC(M_kappa, p_order);

% Store computed PC basis functions
Psi_com = [Psi_com; Psi_u];

%% 4. Smolyak Quadrature for Stochastic Integration
% Determine the number of quadrature points based on PC order
N_integration = ceil((2 * p_order + 1) / 2);

% Generate Smolyak quadrature grid and weights using Gauss-Hermite rule
[xs, ws] = smolyak_grid(M_kappa, N_integration, @gauss_hermite_rule);
xs = xs';   % Ensure dimensions match expectations

%% 5. Generate Stochastic Field Realizations
% Number of stochastic samples from the Smolyak grid
N_xi = size(xs, 1);

% Create vectors for the mean and standard deviation of the field
mean_kappa_vec = mu_kappa * ones(1, nuNodes);
std_kappa_vec = sig_kappa * ones(1, nuNodes);

% Compute realizations of the Gaussian field κ using KL expansion
kappa_i = repmat(mean_kappa_vec, N_xi, 1) + ...
    repmat(std_kappa_vec, N_xi, 1) .* ((eigvec * diag(sqrt(eigval))) * xs')';

% Apply exponential transform for the log-normal Young’s modulus field
E_i = exp(kappa_i)';

%% 6. Assign Results Back to BVP
BVP.ssfem.eigval = eigval;          % KL eigenvalues
BVP.ssfem.eigvec = eigvec;          % KL eigenvectors
BVP.ssfem.M_kappa = M_kappa;        % Number of KL terms
BVP.ssfem.Psi_u = Psi_u;            % Polynomial Chaos basis functions
BVP.ssfem.PsiSqNorm = PsiSqNorm;    % Norms of PC basis functions
BVP.ssfem.P_u = P_u;                % Number of PC terms
BVP.ssfem.Psi_com = Psi_com;        % Combined PC basis functions
BVP.ssfem.xs = xs;                  % Smolyak quadrature points
BVP.ssfem.ws = ws;                  % Smolyak quadrature weights
BVP.ssfem.kappa_i = kappa_i;        % Realizations of κ field
BVP.ssfem.E_i = E_i;                % Realizations of Young’s modulus
BVP.ssfem.N_xi = N_xi;              % Number of stochastic samples

disp('4. KL Expansion of Young modulus is completed.');
end
