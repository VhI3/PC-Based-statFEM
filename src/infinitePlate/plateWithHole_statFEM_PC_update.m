%% plateWithHole_statFEM_PC_update.m
% -------------------------------------------------------------------------
% Statistical FEM: Polynomial Chaos (PC) Update using Kalman Filter
%
% Description:
% This function updates the displacement field and its uncertainty
% using the identified hyperparameters (rho and eta) through the Kalman filter.
% It processes displacements in both x and y directions.
%
% Inputs:
%   BVP : Structure containing:
%         - Covariance matrices (displacement, noise, model discrepancy)
%         - Extended PC coefficients
%         - Observed data
%
% Outputs:
%   BVP : Updated structure with:
%         - Posterior mean displacement fields
%         - Updated displacement covariance matrices
%         - Updated PC coefficients
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0
% -------------------------------------------------------------------------

function BVP = plateWithHole_statFEM_PC_update(BVP)
%% Assign from BVP
% Extracting necessary variables from the BVP structure
C_d_id = BVP.statfem.C_d_id; % Displacement covariance matrix
Hx_full = BVP.obs.Hx_full;
Hy_full = BVP.obs.Hy_full;

rho_id = BVP.statfem.rho_id;
C_e_PC = BVP.statfem.C_e_PC;
C_uu = BVP.ssfem.C_uu;
GDOFs = BVP.preProc.msh.GDOFs; % Global degrees of freedom
nRed = BVP.statfem.nRed; % Number of realizations
u_ext_pc_coef = BVP.statfem.u_ext_pc_cof; % Extended PC coefficients for displacement
e_ext_pc_coef = BVP.statfem.e_ext_pc_coef; % Extended PC coefficients for noise
d_ext_pc_coef = BVP.statfem.d_ext_pc_coef; % Extended PC coefficients for displacement field
P_ext = BVP.statfem.P_ext; % Extended order of PC
PsiSqNorm_ext = BVP.statfem.PsiSqNorm_ext; % Extended Psi square norms
sum_yObs_pc = BVP.statfem.sum_yObs_pc; % Sum of observed data in PC

%% Calculation
H_full = {};
H_full{1} =  Hx_full;
H_full{2} =  Hy_full;

% Compute the combined covariance matrix
sigma_matrix = {};
sigma_matrix{1} = H_full{1} * C_uu{1} * H_full{1}' + (C_d_id{1} + C_e_PC) / (nRed * rho_id{1}^2);
sigma_matrix{2} = H_full{2} * C_uu{2} * H_full{2}' + (C_d_id{2} + C_e_PC) / (nRed * rho_id{2}^2);

% Invert the sigma_matrix using a stable method
inv_sigma_matrix{1} = inv2(sigma_matrix{1});
inv_sigma_matrix{2} = inv2(sigma_matrix{2});

% Compute the Kalman gain
kalman_Gain{1} = C_uu{1} * H_full{1}' * inv_sigma_matrix{1};
kalman_Gain{2} = C_uu{2} * H_full{2}' * inv_sigma_matrix{2};

% Update the extended PC coefficients for displacement using the Kalman gain
de_ext_pc_coef = {};
de_ext_pc_coef{1} = (d_ext_pc_coef{1} + e_ext_pc_coef)/sqrt(nRed);
de_ext_pc_coef{2} = (d_ext_pc_coef{2} + e_ext_pc_coef)/sqrt(nRed);


u_y_ext_pc_coef = {};
u_y_ext_pc_coef{1} = u_ext_pc_coef{1} + (kalman_Gain{1}/rho_id{1}) * (sum_yObs_pc{1} ./ nRed - H_full{1} * rho_id{1}* u_ext_pc_coef{1} - de_ext_pc_coef{1}  );
u_y_ext_pc_coef{2} = u_ext_pc_coef{2} + (kalman_Gain{2}/rho_id{2}) * (sum_yObs_pc{2} ./ nRed - H_full{2} * rho_id{2}* u_ext_pc_coef{2} - de_ext_pc_coef{2}  );

% Computed the corrected PC coefficients for displacement
u_y_ext_pc_coef_corr = {};
u_y_ext_pc_coef_corr{1} = [rho_id{1} * u_y_ext_pc_coef{1}(:, 1) u_y_ext_pc_coef{1}(:,2:end)];
u_y_ext_pc_coef_corr{2} = [rho_id{2} * u_y_ext_pc_coef{2}(:, 1) u_y_ext_pc_coef{2}(:,2:end)];


mu_u_y_pc = {};
mu_u_y_pc{1} = u_y_ext_pc_coef_corr{1}(:, 1);
mu_u_y_pc{2} = u_y_ext_pc_coef_corr{2}(:, 1);


mu_ux_y_pc_Surf         = makeSurf(BVP.preProc.msh.elementNodes,mu_u_y_pc{1});
mu_uy_y_pc_Surf         = makeSurf(BVP.preProc.msh.elementNodes,mu_u_y_pc{2});


% Compute the updated covariance matrix for the displacement
C_uu_y_pc = {};
C_uu_y_pc{1} = zeros(GDOFs/2);
C_uu_y_pc{2} = zeros(GDOFs/2);
for j = 2:P_ext
  C_uu_y_pc{1} = C_uu_y_pc{1} + PsiSqNorm_ext(j) * u_y_ext_pc_coef_corr{1}(:, j) * u_y_ext_pc_coef_corr{1}(:, j)';
  C_uu_y_pc{2} = C_uu_y_pc{2} + PsiSqNorm_ext(j) * u_y_ext_pc_coef_corr{2}(:, j) * u_y_ext_pc_coef_corr{2}(:, j)';
end

std_ux_y_pc       = sqrt(diag( C_uu_y_pc{1}));
std_ux_y_pc_Surf         = makeSurf(BVP.preProc.msh.elementNodes,std_ux_y_pc);

%% Assign back to BVP
BVP.statfem.mu_u_y_pc = mu_u_y_pc; % Updated mean displacement
BVP.statfem.C_uu_y_pc = C_uu_y_pc; % Updated displacement covariance
BVP.statfem.u_y_ext_pc_coef = u_y_ext_pc_coef_corr; % Updated extended PC coefficients
BVP.statfem.mu_ux_y_pc_Surf = mu_ux_y_pc_Surf;
BVP.statfem.mu_uy_y_pc_Surf = mu_uy_y_pc_Surf;
BVP.statfem.std_ux_y_pc_Surf = std_ux_y_pc_Surf;

disp('11. Statistical FEM PC update completed.');
end
