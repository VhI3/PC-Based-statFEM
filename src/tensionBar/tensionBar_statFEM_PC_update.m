function BVP = tensionBar_statFEM_PC_update(BVP)
% TENSIONBAR_STATFEM_PC_UPDATE - Performs the Bayesian update for the displacement field using the Statistical FEM framework.
% This function applies a Gauss-Markov-Kalman filter (GMKF) to update the displacement field
% by assimilating observed data with prior model predictions, incorporating uncertainties.
%
% Inputs (via BVP structure):
%   - C_d_id: Covariance matrix of model-reality mismatch.
%   - C_uu_pc: Prior covariance matrix of displacement.
%   - C_e_PC: Covariance of observation noise.
%   - H, H_full: Projection matrices mapping global displacements to sensor locations.
%   - u_ext_pc_coef, e_ext_pc_coef, d_ext_pc_coef: Polynomial Chaos (PC) coefficients for displacement, noise, and mismatch.
%   - PsiSqNorm_ext: Square norms of extended PC basis functions.
%   - sum_yObs_pc: Sum of observed data projected into PC space.
%   - rho_id: Estimated scaling factor from hyperparameter estimation.
%   - GDof: Global degrees of freedom in the FEM model.
%   - nRed: Number of realizations (data samples).
%   - P_ext: Number of extended PC terms.
%
% Outputs (updated BVP fields):
%   - mu_u_y_pc: Updated mean displacement field.
%   - C_uu_y_pc: Updated covariance of the displacement field.
%   - CI_u_y_pc: 95% confidence interval of updated displacement.
%   - u_y_ext_pc_coef: Updated PC coefficients for the displacement.
%   - pdf_u_y_tip: Probability density function (PDF) of displacement at the tip node.
%   - mu_z_pc: Updated posterior mean at sensor locations.
%   - u_y_MC: Monte Carlo samples of the updated displacement field.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Assign from BVP
% Extract required variables from the BVP structure
C_d_id = BVP.statfem.C_d_id;                   % Model discrepancy covariance matrix
P_u = BVP.ssfem.P_u;                           % Order of Polynomial Chaos for displacement
H_full = BVP.statfem.H_full;                   % Full projection matrix (maps displacement to sensor measurements)
C_uu_pc = BVP.ssfem.C_uu_pc;                   % Prior displacement covariance (from PC expansion)
C_e_PC = BVP.statfem.C_e_PC;                   % Measurement noise covariance (from PC expansion)
GDof = BVP.preProc.GDof;                       % Global degrees of freedom (DOFs)
nRed = BVP.statfem.nRed;                       % Number of observation realizations
rho_id = BVP.statfem.rho_id;                   % Identified scaling factor (hyperparameter)
u_ext_pc_coef = BVP.statfem.u_ext_pc_coef;     % Extended PC coefficients for displacement
e_ext_pc_coef = BVP.statfem.e_ext_pc_coef;     % Extended PC coefficients for noise
d_ext_pc_coef = BVP.statfem.d_ext_pc_coef;     % Extended PC coefficients for model discrepancy
PsiSqNorm_ext = BVP.statfem.PsiSqNorm_ext;     % Square norms of the extended PC basis
sum_yObs_pc = BVP.statfem.sum_yObs_pc;         % Sum of observed data in PC form
P_ext = BVP.statfem.P_ext;                     % Number of extended PC terms
psi_j_ext = BVP.statfem.psi_j_ext;             % Evaluated PC basis functions at sampled points
nMC = BVP.preProc.nMC;                         % Number of Monte Carlo samples

%% Calculation
% Step 1: Compute the total covariance of observed data (including model discrepancy and noise)
sigma_matrix = H_full * C_uu_pc * H_full' + (C_d_id + C_e_PC) / (nRed * rho_id^2);

% Invert the covariance matrix using a stable method
inv_sigma_matrix = inv2(sigma_matrix);  

% Step 2: Kalman Gain calculation
kalman_Gain = C_uu_pc * H_full' * inv_sigma_matrix;

% Step 3: Update the PC coefficients
% Average model discrepancy and noise coefficients
de_ext_pc_coef = (d_ext_pc_coef + e_ext_pc_coef) / sqrt(nRed);  

% Perform the Bayesian update of PC coefficients using Kalman filter formula
u_y_ext_pc_coef = u_ext_pc_coef + (kalman_Gain / rho_id) * ...
    (sum_yObs_pc ./ nRed - H_full * rho_id * u_ext_pc_coef - de_ext_pc_coef);

% Apply correction to ensure consistency with scaling factor (rho_id)
u_y_ext_pc_coef_corr = [rho_id * u_y_ext_pc_coef(:, 1), u_y_ext_pc_coef(:, 2:end)];

% Extract posterior mean displacement (corresponding to the first PC coefficient)
mu_u_y_pc = u_y_ext_pc_coef_corr(:, 1);  

%% Step 4: Covariance and Confidence Interval Computation
% Updated displacement covariance matrix using the corrected PC coefficients
C_uu_y_pc = zeros(GDof);
for j = 2:P_ext
    C_uu_y_pc = C_uu_y_pc + PsiSqNorm_ext(j) * (u_y_ext_pc_coef_corr(:, j) * u_y_ext_pc_coef_corr(:, j)');
end

% Compute the 95% confidence interval
CI_u_y_pc = 1.96 * sqrt(diag(C_uu_y_pc));  

% Alternative covariance computation for verification (comparison with Kalman equations)
C_u_y1 = C_uu_pc - C_uu_pc * H_full' * inv2((C_d_id + C_e_PC)/(rho_id^2 * nRed) + H_full * C_uu_pc * H_full') * H_full * C_uu_pc;
C_u_y3 = (eye(size(C_uu_pc)) - kalman_Gain * H_full) * C_uu_pc;

% Debug information to verify covariance consistency
fprintf("C_u_y1 = %f\n", C_u_y1(end, end));
fprintf("C_u_y3 = %f\n", C_u_y3(end, end));
fprintf("C_u_yp = %f\n", C_uu_y_pc(end, end));

% The posterior mean displacement on sensors
mu_z_pc = H_full*mu_u_y_pc;
C_z_pc = rho_id^2*H_full*C_uu_y_pc*H_full' + C_d_id/sqrt(nRed);
CI_z_pc = sqrt(diag(C_z_pc)) * 1.96;

%% Step 5: Posterior Evaluation at Monte Carlo Samples
% Reconstruct posterior displacement samples using the updated PC coefficients
u_y_MC = zeros(GDof, nMC);  
for j = 1:P_ext
    fprintf('Evaluating PC expansion %g/%g\n', j, P_ext);
    psi_j = psi_j_ext(:, j);      % PC basis evaluated at sample points
    d_j = u_y_ext_pc_coef_corr(:, j);  % Corresponding PC coefficient
    u_y_MC = u_y_MC + d_j * psi_j';    % Accumulate displacement samples
end

% Evaluate displacement at the bar tip and estimate PDF using KDE
u_y_MC_tip = u_y_MC(end, :);                       
[pdf_u_y_tip, x4] = ksdensity(u_y_MC_tip);        

%% Assign Back to BVP
BVP.statfem.mu_u_y_pc = mu_u_y_pc;                  % Posterior mean displacement
BVP.statfem.C_uu_y_pc = C_uu_y_pc;                  % Posterior covariance matrix
BVP.statfem.CI_u_y_pc = CI_u_y_pc;                  % Posterior confidence interval
BVP.statfem.u_y_ext_pc_coef = u_y_ext_pc_coef_corr; % Updated PC coefficients
BVP.statfem.pdf_u_y_tip.pdf_u_y_tip = pdf_u_y_tip;  % Posterior PDF at the bar tip
BVP.statfem.pdf_u_y_tip.u = x4;                     % Support points for the posterior PDF
BVP.statfem.u_y_MC = u_y_MC;                        % Posterior Monte Carlo displacement samples
BVP.statfem.mu_z_pc = mu_z_pc;
BVP.statfem.CI_z_pc = CI_z_pc;
end