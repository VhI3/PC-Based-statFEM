function BVP = tensionBar_statFEM_modelReality_KL_Expansion_converge(BVP, nRed, ld_case, calcase)
% This function performs the Karhunen-Loève (KL) expansion for modeling the 
% model-reality mismatch in the tension bar problem. It uses Gaussian random 
% fields with predefined correlation lengths to represent discrepancies between 
% the computational model and reality. The function reduces the dimension using 
% KL expansion and applies Hermite Polynomial Chaos (PC) expansion.
%
% Inputs:
%   BVP     - Boundary Value Problem structure containing problem settings.
%   nRed    - Number of random samples for reduced order modeling.
%   ld_case - Index for selecting correlation length from predefined vector.
%
% Outputs:
%   Updated BVP structure with sensor eigenvalues, eigenvectors, sampled fields, 
%   and Polynomial Chaos expansions.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Extract Preprocessing Variables from BVP
P_d = BVP.preProc.P_d; % Order of Hermite Polynomial Chaos
L = BVP.preProc.geometry.L; % Length of the domain
tmp_eta_matrix_preDefined = BVP.preProc.tmp_eta_matrix_preDefined; % Predefined eta matrix (KL coefficients)
Psi_com = BVP.statfem.Psi_com; % Composite Psi matrix from statFEM
nMC = BVP.preProc.nMC; % Number of Monte Carlo simulations
psi_j_ext = BVP.ssfem.psi_j_ext; % Extended Psi evaluations

%% Sensor Configuration and Covariance Definition
% Sensor locations and visualization settings
nSen0 = 97; % Number of initial sensors
senCoordinates0 = linspace(2, L - 2, nSen0)'; % Evenly spaced sensors (excluding domain ends)

dotSize = 0.5; % Plotting parameter for sensor size
allp = 0.1; % Transparency for plotting

% Define correlation length for the discrepancy field
corrLength_d_vector = [10, 25, 50, 100];
corrLength_d = corrLength_d_vector(ld_case); % Select correlation length

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
        cov_func_d = @(x1, x2) exp(-0.5 * ((abs(x1 - x2)) / corrLength_d)^2);

    case 'matern_5_2'
        % Matérn Kernel (ν = 5/2):
        % cov(x1, x2) = (1 + sqrt(5)*r + 5r²/3) * exp(-sqrt(5)*r), where r = |x1 - x2|/corrLength_E
        cov_func_d = @(x1, x2) (1 + sqrt(5)*abs(x1 - x2)/corrLength_d + 5/3*(abs(x1 - x2)^2)/(corrLength_d^2)) ...
            .* exp(-sqrt(5)*abs(x1 - x2)/corrLength_d);

    case 'matern_3_2'
        % Matérn Kernel (ν = 3/2):
        % cov(x1, x2) = (1 + sqrt(3)*r) * exp(-sqrt(3)*r)
        cov_func_d = @(x1, x2) (1 + sqrt(3)*abs(x1 - x2)/corrLength_d) ...
            .* exp(-sqrt(3)*abs(x1 - x2)/corrLength_d);

    case 'exponential'
        % Exponential Kernel:
        % cov(x1, x2) = exp(-|x1 - x2| / corrLength_E)
        cov_func_d = @(x1, x2) exp(-abs(x1 - x2)/corrLength_d);

    otherwise
        error('Unknown kernel type selected.');
end


%% Karhunen-Loève Expansion for the Discrepancy Field
% Perform KL expansion using Fredholm-Galerkin approach at full sensor resolution
[sensor_eigVal0, ~, sensor_eigenVec0, sensor_MLterms0] = ...
  KL_Sensor_Galerkin_local_d(senCoordinates0, nSen0, 5, {[senCoordinates0(1), senCoordinates0(end)]}, cov_func_d, 1);

M_d = sensor_MLterms0; % Number of retained KL terms

% Reduce number of sensors for observations
% nSen = 9; % Number of selected sensors for observation
% nSen = 17; % Number of selected sensors for observation
nSen = 33; % Number of selected sensors for observation
senCoordinates_obs = linspace(2, L - 2, nSen); % Reduced sensor coordinates
[~, idxInA] = ismember(senCoordinates_obs, senCoordinates0); % Indices of selected sensors

% Extract sensor information after reduction
senCoordinates = senCoordinates0(idxInA);
sensor_eigVal = sensor_eigVal0;
sensor_eigenVec = sensor_eigenVec0(idxInA, :);

%% Variance Error Analysis of KL Expansion
% Compute variance error of the discretized random field
d_error = ones(nSen, 1) - sum(sensor_eigenVec.^2 * diag(sensor_eigVal), 2);
varianceError = abs(d_error);
meanVarianceError = mean(varianceError); % Average variance error

% Compute eigenvalue-weighted eigenvectors for the discrepancy field
sensor_eigMat = sensor_eigenVec .* sqrt(sensor_eigVal');

%% Polynomial Chaos Expansion for Discrepancy Field
% Extract KL coefficients from predefined eta matrix
eta_matrix_preDefined = tmp_eta_matrix_preDefined(1:M_d);

% Generate Hermite Polynomial Chaos basis functions
[Psi_d, ~] = Hermite_PC_d(M_d, P_d);

% Initialize random samples for discrepancy field
rng(200, 'twister'); % For reproducibility
chi = randn(nRed, M_d); % Reduced random samples
chi_sample = randn(nMC, M_d); % Monte Carlo samples

% Evaluate Hermite PC basis functions at sample points
psi_d = cell(P_d, M_d); % Store evaluated Hermite polynomials
for k = 1:M_d
  chi_n = chi(:, k); % Reduced samples for the k-th KL mode
  chi_v = sym(sprintf('chi_%d', k)); % Symbolic variable for the k-th mode
  chi_sample_n = chi_sample(:, k); % Monte Carlo samples for the k-th mode
  for l = 1:P_d
    psi = Psi_d{l, k}; % Hermite polynomial for mode (k) and order (l)
    psi_d{l, k} = double(subs(psi, chi_v, chi_n)); % Evaluate at reduced samples
    psi_j_ext = [psi_j_ext, double(subs(psi, chi_v, chi_sample_n))]; % Append evaluations at MC samples
  end
end

%% Generate KL-Sampled Discrepancy Field
% Compute the discrepancy field using KL and PC expansions
d_kl_sampled = zeros(nSen, nRed); % Initialize discrepancy field samples
for s = 1:nSen
  for k = 1:M_d
    for l = 1:P_d
      d_kl_sampled(s, :) = d_kl_sampled(s, :) + ...
        sensor_eigMat(s, k) * eta_matrix_preDefined(k, l) * psi_d{l, k}';
    end
  end
end

%% Compute Square Norm Matrix for Hermite PC Basis
ChiSqNorm_matrix = zeros(P_d * M_d); % Initialize norm matrix
counter = 1;
for iM_d = 1:M_d
  for iP_d = 1:P_d
    ChiSqNorm_matrix(counter, counter) = factorial(iP_d); % Norm from factorial of PC order
    counter = counter + 1;
  end
end
ChiSqNorm = diag(ChiSqNorm_matrix); % Extract diagonal norms

% Update composite Psi matrix with reshaped basis functions
Psi_com = [Psi_com; reshape(Psi_d, M_d * P_d, 1)];

%% Assign Results Back to BVP Structure
BVP.statfem.nSen = nSen; % Number of sensors
BVP.statfem.nRed = nRed; % Number of reduced samples
BVP.statfem.senCoordinates = senCoordinates; % Final sensor coordinates
BVP.statfem.sensor_eigVal = sensor_eigVal; % KL eigenvalues
BVP.statfem.M_d = M_d; % Number of KL modes
BVP.statfem.sensor_eigenVec = sensor_eigenVec; % KL eigenvectors
BVP.statfem.sensor_eigMat = sensor_eigMat; % Weighted eigenvectors
BVP.statfem.Psi_d = Psi_d; % Hermite PC basis functions
BVP.statfem.ChiSqNorm_matrix = ChiSqNorm_matrix; % Norm matrix of PC basis
BVP.statfem.ChiSqNorm = ChiSqNorm; % Diagonal norms
BVP.statfem.d_kl_sampled = d_kl_sampled; % Sampled discrepancy field
BVP.statfem.eta_matrix_preDefined = eta_matrix_preDefined; % KL coefficients
BVP.statfem.Psi_com = Psi_com; % Updated composite Psi
BVP.plotSettings.dotSize = dotSize; % Plot dot size
BVP.plotSettings.allp = allp; % Plot transparency
BVP.statfem.chi = chi; % Reduced random variables
BVP.statfem.psi_j_ext = psi_j_ext; % Extended Psi evaluations
BVP.statfem.corrLength_d_vector = corrLength_d_vector; % Correlation length vector
BVP.statfem.varianceError = varianceError; % Variance error from KL expansion
BVP.statfem.corrLength_d = corrLength_d; % Selected correlation length

end
