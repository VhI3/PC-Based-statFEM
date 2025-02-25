function BVP = tensionBar_FEM_preprocess_KL_Expansion(BVP)
% tensionBar_FEM_preprocess_KL_Expansion
% --------------------------------------------------------------
% Performs the Karhunen-Loève (KL) expansion
% for the stochastic modeling of the Young's modulus field in a 1D tension bar problem.
%
% This function:
%   - Computes the KL expansion of a weakly-stationary Gaussian random field representing Young's modulus.
%   - Uses the Fredholm-Galerkin method to calculate eigenvalues and eigenvectors.
%   - Generates stochastic samples using the Smolyak quadrature rule.
%   - Estimates the variance error for the discretized random field.
%
% Inputs:
%   BVP - Boundary Value Problem structure with preprocessing parameters.
%
% Outputs:
%   BVP - Updated structure with:
%         - KL expansion terms (eigenvalues, eigenvectors)
%         - Polynomial chaos basis functions and norms
%         - Stochastic realizations of the field
%         - Variance error of the random field approximation
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Extract parameters from BVP structure
cov_func = BVP.preProc.cov_func;          % Covariance function for the random field
p_order = BVP.preProc.p_order;            % Polynomial chaos expansion order
mu_kappa = BVP.preProc.mu_kappa;          % Mean of the log-random field
sig_kappa = BVP.preProc.sig_kappa;        % Standard deviation of the log-random field
preProc_Variables = BVP.preProc;          % Preprocessing variables
GDof = preProc_Variables.GDof;            % Global degrees of freedom

% Initialize Psi combination matrix
Psi_com = [];

%% KL Expansion using Fredholm-Galerkin method
% Computes the eigenvalues (eigval) and eigenvectors (eigvec) for the random field.
[eigVal_numerics, ~, eigenVec_numerics, M_kappa_numerics] = KL_fredholm_Galerkin_local_1D(preProc_Variables, cov_func);

eigval = eigVal_numerics;     % KL eigenvalues
eigvec = eigenVec_numerics;   % KL eigenvectors
M_kappa = M_kappa_numerics;   % Number of KL terms retained

%% Variance error calculation
% Checks the approximation accuracy of the discretized random field.
varianceError = abs(ones(GDof, 1) - sum((eigvec.^2) * diag(eigval), 2));
meanVarianceError = mean(varianceError);

%% Polynomial Chaos Expansion (PCE) calculations
% Generate polynomial chaos basis functions and their norms.
[~, Psi_u, ~, PsiSqNorm, P_u] = Hermite_PC(M_kappa, p_order);

% Store Psi basis functions for combined use in later calculations
Psi_com = [Psi_com; Psi_u];

%% Pseudo Spectral Method for Stochastic Expansion
% Generate Smolyak quadrature points and weights for integration.
N_integration = ceil((2 * p_order + 1) / 2);  % Number of integration points
[xs, ws] = smolyak_grid(M_kappa, N_integration, @gauss_hermite_rule);
xs = xs';                                     % Ensure correct dimension

% Generate stochastic realizations of the log-random field kappa
N_xi = size(xs, 1);                           % Number of stochastic samples
mean_kappa_vec = mu_kappa * ones(1, GDof);    % Mean kappa vector
std_kappa_vec = sig_kappa * ones(1, GDof);    % Standard deviation vector

% Stochastic samples of kappa using KL expansion and PCE
kappa_i = repmat(mean_kappa_vec, N_xi, 1) + ...
    repmat(std_kappa_vec, N_xi, 1) .* ((eigvec * diag(sqrt(eigval))) * xs')';

E_i = exp(kappa_i);                           % Young’s modulus realizations

%% Assign results back to BVP structure
BVP.ssfem.M_kappa = M_kappa;                  % Number of KL terms
BVP.ssfem.eigval = eigval;                    % KL eigenvalues
BVP.ssfem.eigvec = eigvec;                    % KL eigenvectors
BVP.ssfem.P_u = P_u;                          % PCE basis functions
BVP.ssfem.PsiSqNorm = PsiSqNorm;              % Norm of PCE basis functions
BVP.ssfem.Psi_u = Psi_u;                      % PCE polynomials
BVP.ssfem.N_integration = N_integration;      % Number of integration points
BVP.ssfem.xs = xs;                            % Smolyak grid points
BVP.ssfem.ws = ws;                            % Smolyak grid weights
BVP.ssfem.N_xi = N_xi;                        % Number of stochastic samples
BVP.ssfem.mean_kappa_vec = mean_kappa_vec;    % Mean of kappa
BVP.ssfem.std_kappa_vec = std_kappa_vec;      % Standard deviation of kappa
BVP.ssfem.kappa_i = kappa_i;                  % Stochastic samples of kappa
BVP.ssfem.E_i = E_i;                          % Exponential field (Young’s modulus)
BVP.ssfem.varianceError = varianceError;      % Variance error
BVP.ssfem.meanVarianceError = meanVarianceError; % Mean variance error

% Store combined Psi basis
BVP.statfem.Psi_com = Psi_com;

disp('3. KL expansion of Young modulus is completed.');
end
