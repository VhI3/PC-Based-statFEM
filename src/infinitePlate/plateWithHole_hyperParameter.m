%% plateWithHole_hyperParameter.m
% -------------------------------------------------------------------------
% Statistical FEM: Hyperparameter Estimation for Model-Reality Mismatch
%
% Description:
% This function estimates hyperparameters (scaling factor `rho` and KL expansion
% coefficients `eta`) to align the computational model with observed data.
% The estimation uses a negative log-likelihood function and optimizes parameters
% for both x- and y-direction sensor observations.
%
% Inputs:
%   BVP : Boundary Value Problem structure containing:
%         - Observed data
%         - Stochastic displacement samples
%         - Polynomial Chaos and KL expansions
%
% Outputs:
%   BVP : Updated structure containing:
%         - Identified scaling factors (rho)
%         - Identified KL expansion coefficients (eta)
%         - Updated covariance matrices
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0
% -------------------------------------------------------------------------

function BVP = plateWithHole_hyperParameter(BVP)
%% Assign from BVP
% Extracting necessary variables from the BVP structure for hyperparameter estimation
% eta_matrix_preDefined = BVP.statfem.eta_matrix_preDefined; % Predefined eta matrix
% rho_preDefined = BVP.preProc.rho_preDefined; % Predefined scaling factor rho
N_xi = BVP.ssfem.N_xi; % Number of stochastic samples for displacement
u_xii = BVP.ssfem.u_xii; % Displacement samples
ws = BVP.ssfem.ws; % Weights for stochastic samples
% H = BVP.obs.H_full; % Projection matrix
C_e_PC = BVP.statfem.C_e_PC; % Polynomial Chaos expansion of noise covariance
ChiSqNorm_matrix = BVP.statfem.ChiSqNorm_matrix; % Square norm matrix for Chi
sensor_eigMat = BVP.statfem.sensor_eigMat; % Eigenvalue matrix for sensor
nSen = BVP.preProc.statfem.nuSensors; % Number of sensors
nRed = BVP.statfem.nRed; % Number of realizations
yObs = BVP.statfem.Y_exp; % Observed data
P_u = BVP.ssfem.P_u; % Order of PC expansion for u
M_d = BVP.statfem.M_d; % Number of terms in KL expansion for d
P_d = 1; % Order of PC expansion for d
P_ext = BVP.statfem.P_ext; % Extended order of PC
P = BVP.obs.P;

%% Calculation


options = optimoptions('fminunc', ...
    'Display', 'iter',...
    'Algorithm', 'quasi-newton', ...
    'FiniteDifferenceType', 'central', ...
    'FunctionTolerance', 1e-3, ...
    'StepTolerance', 1e-3, ...
    'MaxIterations', 100);




% Starting point for optimization includes initial rho and log of predefined eta values
% modelReality_hyperparameter_startPoint = log(eta_matrix_preDefined');
rng(10,"twister")
modelReality_hyperparameter_startPoint = log(1+rand(M_d,1))';
startPar = [1, modelReality_hyperparameter_startPoint];

hyperparameters_id_2D = zeros(2,M_d+1);
for ii = 1:2
    % Negative log-likelihood function for observed data
    neg_loglike_Yobs = @(hyperparameters_id) ...
        neg_loglikelihood_Yobs2D(hyperparameters_id, N_xi, nRed, nSen, C_e_PC, yObs(ii:2:end,:), P{ii}, u_xii{ii}, sensor_eigMat, ChiSqNorm_matrix, ws);

    % Minimize the negative log-likelihood to estimate hyperparameters
    hyperparameters_id_2D(ii,:) = fminunc(neg_loglike_Yobs, startPar, options);
end
hyperparameters_id = {};
hyperparameters_id{1} = zeros(2,M_d+1);
hyperparameters_id{2} = zeros(2,M_d+1);
for ii = 1:2
    hyperparameters_id{ii} = [hyperparameters_id_2D(ii,1) sort(exp(hyperparameters_id_2D(ii,2:end)),'descend')];
end



% Extract the identified rho and eta matrix from the optimized hyperparameters
rho_id = {};
rho_id{1} = hyperparameters_id{1}(1);
rho_id{2} = hyperparameters_id{2}(1);
eta_matrix_klweise_id = {};
eta_matrix_klweise_id{1} = diag(hyperparameters_id{1}(2:end));
eta_matrix_klweise_id{2} = diag(hyperparameters_id{2}(2:end));

% Compute the extended displacement coefficients using the identified eta matrix
d_ext_pc_coef = {};
d_ext_pc_coef{1} = zeros(nSen, P_ext);
d_ext_pc_coef{2} = zeros(nSen, P_ext);

d_ext_pc_coef{1}(:, P_u + 1:P_u + M_d * P_d) = sensor_eigMat * eta_matrix_klweise_id{1};
d_ext_pc_coef{2}(:, P_u + 1:P_u + M_d * P_d) = sensor_eigMat * eta_matrix_klweise_id{2};


% Calculate the covariance matrix for the identified displacement field
C_d_id = {};
C_d_id{1} = (sensor_eigMat * eta_matrix_klweise_id{1}) * ChiSqNorm_matrix * (eta_matrix_klweise_id{1}' * sensor_eigMat');
C_d_id{2} = (sensor_eigMat * eta_matrix_klweise_id{2}) * ChiSqNorm_matrix * (eta_matrix_klweise_id{2}' * sensor_eigMat');

%
%% Assign back to BVP
BVP.statfem.rho_id = rho_id;
BVP.statfem.hyperparameters_id = hyperparameters_id;
BVP.statfem.d_ext_pc_coef = d_ext_pc_coef;
BVP.statfem.C_d_id = C_d_id;
% BVP.statfem.rhoFeed = rhoFeed;
end

function [f] = neg_loglikelihood_Yobs2D(params, N_xi, nRed, nSen, C_e_PC, y, H, u, sensor_eigMat, ChiSqNorm_matrix, w)
% This function calculates the negative log-likelihood of observed data given the model parameters.
% It is used to estimate the model parameters by minimizing this value.

% Inputs:
% params: vector of parameters (rho and ln_sig_d)
% N_xi: Number of stochastic samples for displacement
% nRed: Number of realizations
% nSen: Number of sensors
% C_e_PC: Polynomial Chaos expansion of noise covariance
% y: Observed data
% H: Projection matrix
% u: Displacement samples
% sensor_eigMat: Eigenvalue matrix for sensor
% ChiSqNorm_matrix: Square norm matrix for Chi
% w: Weights for stochastic samples

% Output:
% f: Negative log-likelihood value

% Extract rho and ln_sig_d from the parameters
rho = params(1);
ln_sig_d = params(2:end);

% Construct the eta matrix from ln_sig_d
eta_matrix_klweise = diag(exp(ln_sig_d'));

% Compute the covariance matrix for the displacement field using the eta matrix
C_d_PC = (sensor_eigMat * eta_matrix_klweise) * ChiSqNorm_matrix * (eta_matrix_klweise' * sensor_eigMat');
% C_d_PC = [C_d_PC0 zeros(size(C_d_PC0,1));
%    zeros(size(C_d_PC0,1)) C_d_PC0];
% Sum the displacement and noise covariance matrices
C_d_e = C_d_PC + C_e_PC;

% % % Perform Cholesky decomposition for numerical stability
L = chol(C_d_e, 'lower');
% Solve for the inverse of C_d_e using Cholesky factor
% inv_C_d_e = chol_solve(L, eye(nSen));
L_inv = L \ eye(size(L));
inv_C_d_e = L_inv' * L_inv;
% Calculate determinant using Cholesky factor
tmp_det = 2 * nRed * sum(log(diag(L)));

% inv_C_d_e = inv(C_d_e);
% tmp_det = det(C_d_e);


% Initialize the log-likelihood function
% f1 = -0.5 * tmp_det - 0.5 * nRed * nSen * log(2 * pi);
f1 = 0.5 * tmp_det + 0.5 * nRed * nSen * log(2 * pi);

% Initialize y_star for accumulating the weighted sum of squared differences
y_star = zeros(N_xi, 1);
for j = 1:N_xi
    for i = 1:nRed
        % y_star(j, 1) = y_star(j, 1) + (-0.5) * (y(:, i) - rho * H * u(:, j))' * (inv_C_d_e) * (y(:, i) - rho * H * u(:, j)) * w(j);
        y_star(j, 1) = y_star(j, 1) + (-0.5) * (y(:, i) - rho * H * u(:, j))' * (inv_C_d_e) * (y(:, i) - rho * H * u(:, j));
    end
end

% Compute the log-sum-exp for numerical stability
% [f2_, ~] = logsumexp(y_star);
[f2_] = logsumexp_weighted(y_star,w);
% Calculate the final negative log-likelihood
f = (f1 - f2_);
end


function [lse,sm] = logsumexp(x)
%LOGSUMEXP  Log-sum-exp function.
%    lse = LOGSUMEXP(x) returns the log-sum-exp function evaluated at
%    the vector x, defined by lse = log(sum(exp(x)).
%    [lse,sm] = LOGSUMEXP(x) also returns the softmax function evaluated
%    at x, defined by sm = exp(x)/sum(exp(x)).
%    The functions are computed in a way that avoids overflow and
%    optimizes numerical stability.


if ~isvector(x), error('Input x must be a vector.'), end

n = length(x);
s = 0; e = zeros(n,1);
[xmax,k] = max(x); a = xmax;
s = 0;
for i = 1:n
    e(i) = exp(x(i)-xmax);
    if i ~= k
        s = s + e(i);
    end
end
lse = a + log1p(s);
if nargout > 1
    sm = e/(1+s);
end

end


function [lse] = logsumexp_weighted(x, w)
%LOGSUMEXP_WEIGHTED Log-sum-exp function with weights.
%    lse = LOGSUMEXP_WEIGHTED(x, w) returns the log-sum-exp function
%    evaluated at the vector x with weights w, defined by
%    lse = log(sum(w .* exp(x))).
%    The function is computed in a way that avoids overflow and
%    optimizes numerical stability.

% Check if the input x is a vector and if w has the same length as x
if ~isvector(x)
    error('Input x must be a vector.');
end
if length(w) ~= length(x)
    error('Input w must have the same length as x.');
end

% Get the length of the vector x
n = length(x);

% Initialize an array e with zeros
e = zeros(n, 1);

% Find the maximum value in x and its index
[xmax, k] = max(x);

% Set the shift value a to the maximum value
a = xmax;

% Initialize the sum s to zero
s = 0;

% Compute the exponentials and the weighted sum
for i = 1:n
    e(i) = w(i) * exp(x(i) - xmax);
    s = s + e(i);
end

% Compute the log-sum-exp value using the log1p function for better accuracy
lse = a + log(s); % changed from log1p to log
end