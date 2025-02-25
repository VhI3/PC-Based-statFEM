%% plateWithHole_SSFEM_process_LE_PC_Expansion.m
% --------------------------------------------------------------
% Stochastic Spectral Finite Element Method (SSFEM) with Polynomial Chaos (PC) Expansion
%
% Description:
% This function computes the stochastic displacement field for a 2D plate with a hole
% under linear elasticity using SSFEM and PC expansion. The method involves:
%   1. Realizing displacement fields for each stochastic sample of Young’s modulus.
%   2. Projecting these realizations onto polynomial chaos basis functions.
%   3. Computing the mean, covariance, and confidence intervals of displacements.
%
% Inputs:
%   BVP : Boundary Value Problem structure with preprocessing variables.
%
% Outputs:
%   BVP : Updated structure containing SSFEM results: PC coefficients, mean field, covariance, and MC samples.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

function BVP = plateWithHole_SSFEM_process_LE_PC_Expansion(BVP)

%% Extracting necessary variables from BVP structure
N_xi                        = BVP.ssfem.N_xi; % Number of stochastic samples
preProc_Variables           = BVP.preProc; % Preprocessing variables
E_vector                    = BVP.ssfem.E_i; % Young's modulus realizations
xs                          = BVP.ssfem.xs; % Smolyak grid points
ws                          = BVP.ssfem.ws; % Smolyak weights
M_kappa                     = BVP.ssfem.M_kappa; % Number of terms in KL expansion
GDOFs                       = BVP.preProc.msh.GDOFs; % Degrees of freedom
P_u                         = BVP.ssfem.P_u; % Number of polynomial chaos terms
Psi_u                       = BVP.ssfem.Psi_u; % Polynomial chaos basis functions
activeDofs                  = BVP.preProc.BC.activeDOFs; % Active degrees of freedom
PsiSqNorm                   = BVP.ssfem.PsiSqNorm; % Norms of the polynomial chaos basis functions
nMC                         = BVP.preProc.nMC;
%
%
%
%% Calculation of displacement realizations
u_pc_realization = zeros(GDOFs, N_xi); % Initialize displacement realizations matrix
for i = 1:N_xi
    % Solve for displacement for each realization of E
    [u_pc_realization(:, i), ~] = FEM_plateWithHole_Solver(preProc_Variables, E_vector(:,i));
end

%% Preparation for Polynomial Chaos expansion
% Create symbolic and numeric variables for Polynomial Chaos
xi_v = cell(1,M_kappa);   % xi symbolic vars
xi_n = cell(1,M_kappa);   % xi numeric
for j = 1:M_kappa
    xi_v{j} = sym(sprintf('xi_%d', j)); % symbolic variables
    xi_n{j} = xs(:, j); % numeric variables (samples)
end

%% Polynomial Chaos coefficients calculation
u_pc_coef = zeros(GDOFs, P_u); % Initialize PC coefficients matrix

for i = 1:P_u
    psi = Psi_u{i}; % Get the i-th polynomial chaos basis function
    % psi_j = double(subs(psi, xi_v, xi_n)); % Evaluate the basis function at the sample points
    
    psi_j_mat(:,i) = double(subs(psi, xi_v, xi_n)); % Evaluate the basis function at the sample points
    psi_j = psi_j_mat(:,i);
    
    % Calculate the coefficients by projection
    for j = 1:N_xi
        u_pc_coef(:, i) = u_pc_coef(:, i) + u_pc_realization(:, j) * psi_j(j) * ws(j);
    end
    
    % Normalize the coefficients
    u_pc_coef(:, i) = u_pc_coef(:, i) / PsiSqNorm(i);
end

%% Extracting results for active degrees of freedom
u_pc_coef_active        = u_pc_coef(activeDofs, :);
mu_u_pc                 = u_pc_coef(:, 1); % Mean displacement
mu_u_pc_active          = mu_u_pc(activeDofs); % Mean displacement for active DOFs

mu_ux_pc                = mu_u_pc(1:2:end,:);
mu_uy_pc                = mu_u_pc(2:2:end,:);
mu_ux_pc_Surf           = makeSurf(BVP.preProc.msh.elementNodes,mu_ux_pc);
mu_uy_pc_Surf           = makeSurf(BVP.preProc.msh.elementNodes,mu_uy_pc);

%% Covariance of the displacement response
C_uu_pc = zeros(GDOFs); % Initialize covariance matrix

% Calculate the covariance matrix by summing over the outer products of the PC coefficients
for j = 2:P_u
    C_uu_pc = C_uu_pc + PsiSqNorm(j) * u_pc_coef(:, j) * u_pc_coef(:, j)';
end

% Symmetrize the covariance matrix
C_uu_pc = 0.5 * (C_uu_pc + C_uu_pc');
CI_u_pc = sqrt(diag(C_uu_pc)) * 1.96; % 95% confidence interval
C_uu_pc_active = C_uu_pc(2:end, 2:end); % Covariance for active DOFs
% cov for x
C_uux_pc        = C_uu_pc(1:2:end,1:2:end);
std_ux_pc       = sqrt(diag(C_uux_pc));
CI_ux_pc        = std_ux_pc*1.96;
% cov for y
C_uuy_pc        = C_uu_pc(2:2:end,2:2:end);
std_uy_pc       = sqrt(diag(C_uuy_pc));
CI_uy_pc        = std_uy_pc*1.96;

C_uu = {};
C_uu{1} = C_uux_pc;
C_uu{2} = C_uuy_pc;

std_ux_pc_Surf         = makeSurf(BVP.preProc.msh.elementNodes, std_ux_pc);
std_uy_pc_Surf         = makeSurf(BVP.preProc.msh.elementNodes, std_uy_pc);


% Reconstruct displacement field for each stochastic sample
u_xi = u_pc_coef * psi_j_mat.';
ux_xi = u_xi(1:2:end,:);
uy_xi = u_xi(2:2:end,:);
u_xii{1} = ux_xi;
u_xii{2} = uy_xi;
uMC  = zeros(GDOFs, nMC);
rng(100, 'twister') % Setting the random number generator for reproducibility
xi = randn(nMC,M_kappa);
for j = 1:M_kappa
    xi_mc{j} = xi(:, j); % numeric variables (samples)
end
for j = 1:P_u
    msg = fprintf('Evaluating PC expansion %g/%g', j, P_u);
    % using syms
    psi   = Psi_u{j};
    psi_j = double(subs(psi, xi_v, xi_mc));
    psi_j_ext(:,j) = psi_j;
    d_j   = u_pc_coef(:, j);
    for i = 1:GDOFs
        uMC(i,:) = uMC(i,:) + d_j(i)*psi_j';
    end
    fprintf(repmat('\b',1,msg));
end


%% Assigning results back to BVP structure
BVP.ssfem.u_pc_coef             = u_pc_coef;
BVP.ssfem.u_pc_coef_active      = u_pc_coef_active;
BVP.ssfem.mu_u_pc               = mu_u_pc;
BVP.ssfem.mu_u_pc_active        = mu_u_pc_active;
BVP.ssfem.mu_ux_pc              = mu_ux_pc;
BVP.ssfem.mu_uy_pc              = mu_uy_pc;
BVP.ssfem.mu_ux_pc_Surf         = mu_ux_pc_Surf;
BVP.ssfem.mu_uy_pc_Surf         = mu_uy_pc_Surf;
BVP.ssfem.C_uu_pc               = C_uu_pc;
BVP.ssfem.C_uu                  = C_uu;
BVP.ssfem.C_uu_pc_active        = C_uu_pc_active;
BVP.ssfem.CI_u_pc               = CI_u_pc;
BVP.ssfem.C_uux_pc              = C_uux_pc;
BVP.ssfem.C_uuy_pc              = C_uuy_pc;
BVP.ssfem.std_ux_pc             = std_ux_pc;
BVP.ssfem.std_uy_pc             = std_uy_pc;
BVP.ssfem.CI_uy_pc              = CI_uy_pc;
BVP.ssfem.CI_ux_pc              = CI_ux_pc;
BVP.ssfem.std_ux_pc_Surf        = std_ux_pc_Surf;
BVP.ssfem.std_uy_pc_Surf        = std_uy_pc_Surf;
BVP.ssfem.u_xi                  = u_xi;
BVP.ssfem.ux_xi = ux_xi;
BVP.ssfem.uy_xi = uy_xi;
BVP.ssfem.u_xii = u_xii;
BVP.statfem.psi_j_ext = psi_j_ext;
BVP.ssfem.uMC = uMC;

disp('5. PC expansion of displacement is completed.');
end