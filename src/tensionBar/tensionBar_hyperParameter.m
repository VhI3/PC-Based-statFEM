function BVP = tensionBar_hyperParameter(BVP)
% tensionBar_hyperParameter - Estimates hyperparameters (rho and eta) using maximum likelihood.
%
% Description:
%   This function estimates the scaling factor (rho) and model-reality mismatch hyperparameters (eta)
%   by minimizing the negative log-likelihood of the observed data given the model predictions.
%   It uses the quasi-Newton optimization algorithm through MATLAB's `fminunc` function.
%
% Inputs:
%   BVP - Boundary Value Problem structure with precomputed data and observations.
%
% Outputs (updated fields in BVP.statfem):
%   rho_id                  - Identified scaling factor rho.
%   hyperparameters_id      - Identified hyperparameters vector (rho and eta terms).
%   d_ext_pc_coef           - Extended Polynomial Chaos coefficients for model-reality mismatch.
%   C_d_id                  - Covariance matrix of the identified displacement field.
%   rhoFeed                 - Initial guess for rho (mean of observed-to-predicted ratio).
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Extract Necessary Variables from BVP
eta_matrix_preDefined  = BVP.statfem.eta_matrix_preDefined;  % Predefined eta values (model mismatch)
N_xi                   = BVP.ssfem.N_xi;                     % Number of stochastic samples for displacement
u_xi                   = BVP.ssfem.u_xi;                     % Displacement samples (from PC expansions)
ws                     = BVP.ssfem.ws;                       % Weights for stochastic quadrature points
H                      = BVP.statfem.H;                      % Projection matrix from global DOFs to sensor space
C_e_PC                 = BVP.statfem.C_e_PC;                 % Covariance matrix of the measurement noise
ChiSqNorm_matrix       = BVP.statfem.ChiSqNorm_matrix;       % Square norm matrix for Chi-based PC basis
sensor_eigMat          = BVP.statfem.sensor_eigMat;          % Eigenfunctions (KL basis) for model-reality mismatch
nSen                   = BVP.statfem.nSen;                   % Number of sensors
nRed                   = BVP.statfem.nRed;                   % Number of realizations
yObs                   = BVP.statfem.yObs;                   % Observed data
P_u                    = BVP.ssfem.P_u;                      % Number of PC terms for displacement
M_d                    = BVP.statfem.M_d;                    % Number of KL terms for model-reality mismatch
P_d                    = BVP.preProc.P_d;                    % Number of PC terms for model-reality mismatch
P_ext                  = BVP.statfem.P_ext;                  % Total number of PC terms in the extended basis
mean_yObs              = BVP.statfem.mean_yObs;              % Mean of observed data
U_analytic_active      = BVP.analyticalSolution.U_analytic_active;  % Analytical solution for displacement (active DOFs)

%% Step 1: Estimate Initial Guess for rho
% Project the analytical displacement solution to sensor locations
projected_analytic_u_active = H * U_analytic_active;

% Calculate initial rho guess as the ratio between mean observed and predicted values
rhoVector = mean_yObs ./ projected_analytic_u_active; 
rhoFeed = mean(rhoVector);  % Initial rho estimate

%% Step 2: Set Initial Parameters for Optimization
% Starting point includes:
%   - rho: scaling factor (initialized to 1 for neutrality).
%   - eta: log-transformed predefined eta values for better optimization stability.
modelReality_hyperparameter_startPoint = log(eta_matrix_preDefined');  % Log-transform eta
startPar = [1, modelReality_hyperparameter_startPoint];                % Combined initial guess

%% Step 3: Define the Negative Log-Likelihood Function
% The objective function evaluates how likely the observed data is given the hyperparameters.
neg_loglike_Yobs = @(hyperparameters_id) ...
    neg_loglikelihood_Yobs(hyperparameters_id, N_xi, nRed, nSen, C_e_PC, ...
                           yObs, H, u_xi(2:end, :), sensor_eigMat, ChiSqNorm_matrix, ws);

%% Step 4: Perform Hyperparameter Optimization
% Use the 'quasi-newton' method via fminunc to minimize the negative log-likelihood.
options = optimoptions(@fminunc, 'Display', 'iter', 'Algorithm', 'quasi-newton');
hyperparameters_id = fminunc(neg_loglike_Yobs, startPar, options);

% Post-process:
%   - First element is rho (scaling factor).
%   - Subsequent elements are eta values (retransformed from log scale and sorted).
hyperparameters_id = [hyperparameters_id(1), sort(exp(hyperparameters_id(2:end)), 'descend')];

%% Step 5: Display the Identified Hyperparameters
fprintf("------------------------------------------------\n")
fprintf("Initial rho = %f \n", rhoFeed);
fprintf("Identified rho = %f \n", hyperparameters_id(1));

fprintf("Initial eta values: ");
fprintf("%f \t", eta_matrix_preDefined);
fprintf("\n");

fprintf("Identified eta values: ");
fprintf("%f \t", hyperparameters_id(2:end));
fprintf("\n------------------------------------------------\n");

%% Step 6: Reconstruct Identified Model-Reality Mismatch
rho_id = hyperparameters_id(1);                                     % Identified rho
eta_matrix_klweise_id = diag(hyperparameters_id(2:end));            % Diagonal eta matrix

% Reconstruct Polynomial Chaos coefficients for the model-reality mismatch displacement field
d_ext_pc_coef = zeros(nSen, P_ext);
d_ext_pc_coef(:, P_u + 1:P_u + M_d * P_d) = sensor_eigMat * eta_matrix_klweise_id;

% Calculate the covariance matrix for the identified model-reality mismatch
C_d_id = sensor_eigMat * eta_matrix_klweise_id * ChiSqNorm_matrix * (eta_matrix_klweise_id' * sensor_eigMat');

%% Step 7: Assign Results Back to BVP Structure
BVP.statfem.rho_id                = rho_id;                               % Identified rho
BVP.statfem.hyperparameters_id    = [rho_id; hyperparameters_id(2:end)']; % Full hyperparameter set
BVP.statfem.d_ext_pc_coef         = d_ext_pc_coef;                        % PC coefficients of mismatch
BVP.statfem.C_d_id                = C_d_id;                               % Covariance of mismatch
BVP.statfem.rhoFeed               = rhoFeed;                              % Initial rho estimate

end




function [f] = neg_loglikelihood_Yobs(params, N_xi, nRed, nSen, C_e_PC, y, H, u, sensor_eigMat, ChiSqNorm_matrix, w)
% neg_loglikelihood_Yobs - Calculates the negative log-likelihood of observed data given model parameters.
%
% Description:
%   This function evaluates the likelihood of observed sensor data under the model, given the hyperparameters:
%   - rho: Scaling factor for the displacement predictions.
%   - eta: Model-reality mismatch hyperparameters represented through KL expansion.
%
%   It is used within an optimization routine (e.g., fminunc) to estimate rho and eta by minimizing this value.
%   The formulation assumes a Gaussian model for the observations with combined covariance from model error and noise.
%
% Inputs:
%   params             - Vector of parameters (rho and ln_sig_d).
%   N_xi               - Number of stochastic samples for displacement field.
%   nRed               - Number of realizations (Monte Carlo samples for observations).
%   nSen               - Number of sensors.
%   C_e_PC             - Noise covariance matrix from Polynomial Chaos expansion.
%   y                  - Observed sensor data matrix (size: nSen x nRed).
%   H                  - Projection matrix (maps global displacements to sensor space).
%   u                  - Displacement field samples (size: active DOFs x N_xi).
%   sensor_eigMat      - Matrix of eigenfunctions (KL basis for model-reality mismatch).
%   ChiSqNorm_matrix   - Square norm matrix for Hermite Polynomial Chaos basis.
%   w                  - Quadrature weights for stochastic integration.
%
% Output:
%   f - Negative log-likelihood value (to be minimized during hyperparameter estimation).
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Step 1: Extract Hyperparameters
rho = params(1);                  % Scaling factor for model predictions
ln_sig_d = params(2:end);         % Logarithm of eta hyperparameters (ensures positivity after exp transform)

%% Step 2: Construct Model-Reality Mismatch Covariance
% eta_matrix_klweise: Diagonal matrix of eta values controlling magnitude of each KL term
eta_matrix_klweise = diag(exp(ln_sig_d')); 

% C_d_PC: Covariance matrix of the model-reality mismatch field
C_d_PC = (sensor_eigMat * eta_matrix_klweise) * ChiSqNorm_matrix * (eta_matrix_klweise' * sensor_eigMat');

%% Step 3: Total Covariance (Model Mismatch + Observation Noise)
% Combine the covariance of the mismatch (C_d_PC) and the measurement noise (C_e_PC)
C_d_e = C_d_PC + C_e_PC;

%% Step 4: Cholesky Decomposition (for Efficient Inversion and Determinant Computation)
% Cholesky factorization: C_d_e = L * L', where L is lower triangular
L = chol(C_d_e, 'lower'); 

% Efficient inverse calculation using Cholesky factors
L_inv = L \ eye(size(L));           % Solve L * X = I
inv_C_d_e = L_inv' * L_inv;         % Inverse of C_d_e using L_inv

% Determinant of C_d_e using Cholesky factors: det(C_d_e) = (prod(diag(L)))^2
tmp_det = 2 * nRed * sum(log(diag(L))); 

%% Step 5: Likelihood Terms
% Log-likelihood constant terms (normalizing factor of multivariate Gaussian)
f1 = 0.5 * tmp_det + 0.5 * nRed * nSen * log(2 * pi); 

%% Step 6: Data Misfit Computation
% y_star accumulates weighted squared differences between observations and model predictions
y_star = zeros(N_xi, 1);  

% Iterate over stochastic displacement samples and realizations to accumulate the misfit
for j = 1:N_xi
    for i = 1:nRed
        residual = y(:, i) - rho * H * u(:, j);     % Observation-model residual
        y_star(j) = y_star(j) + (-0.5) * residual' * inv_C_d_e * residual;  % Quadratic form in Gaussian likelihood
    end
end

%% Step 7: Numerical Stability via Log-Sum-Exp Trick
% Convert the accumulated misfit to log-likelihood using weighted log-sum-exp for stability
f2_ = logsumexp_weighted(y_star, w);  

%% Step 8: Negative Log-Likelihood
% The optimization aims to minimize this value (hence the sign conventions)
f = f1 - f2_;  

end


function [lse] = logsumexp_weighted(x, w)
% LOGSUMEXP_WEIGHTED - Computes the weighted log-sum-exp function in a numerically stable way.
%
% Syntax:
%   lse = logsumexp_weighted(x, w)
%
% Description:
%   This function calculates the value of:
%
%      lse = log(sum(w .* exp(x)))
%
%   directly, but in a way that prevents numerical overflow and improves stability.
%   The standard computation of exp(x) can lead to large exponentials, causing floating-point issues.
%   By factoring out the maximum element of x, this function ensures stability.
%
% Inputs:
%   x - Input vector (real values), typically representing log-probabilities or exponent terms.
%   w - Corresponding weight vector (same length as x), scaling each exponential term.
%
% Output:
%   lse - The weighted log-sum-exp value (scalar).
%
% Example:
%   x = [-1000, -2, -1]; w = [0.2, 0.5, 0.3];
%   result = logsumexp_weighted(x, w); % Stable computation of log(sum(w .* exp(x)))
%
% Notes:
%   - If weights are all ones (w = ones(size(x))), this function reduces to the standard log-sum-exp.
%   - logsumexp is commonly used in statistical models, especially with log-likelihoods and posterior computations.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Input Validation
if ~isvector(x)
    error('Input x must be a vector.');
end

if length(w) ~= length(x)
    error('Input w must have the same length as x.');
end

%% Step 1: Numerical Stability via Max Shifting
% Find the maximum element in x to use as a shift constant (avoids large exponentials)
[xmax, ~] = max(x);  
a = xmax;  % Shift value

%% Step 2: Exponentiate Shifted Values and Compute Weighted Sum
% Initialize an array for exponentials and a sum accumulator
n = length(x);
e = zeros(n, 1);  % Exponentials storage
s = 0;            % Weighted sum accumulator

% Calculate weighted exponentials with the shift applied
for i = 1:n
    e(i) = w(i) * exp(x(i) - a);  % w_i * exp(x_i - xmax)
    s = s + e(i);                 % Accumulate sum
end

%% Step 3: Final Log-Sum-Exp Calculation
% Combine the shift and the weighted sum:
% log(sum(w .* exp(x))) = a + log(sum(w .* exp(x - a)))
lse = a + log(s);  

end
