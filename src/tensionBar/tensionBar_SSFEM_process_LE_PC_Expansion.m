function BVP = tensionBar_SSFEM_process_LE_PC_Expansion(BVP)
% This function processes the linear elastic boundary value problem for a tension bar
% using Stochastic Spectral Finite Element Method (SSFEM) with Polynomial Chaos (PC) Expansion.
% It calculates the mean and covariance of the displacement field and updates the BVP structure.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Extracting necessary variables from BVP structure
N_xi = BVP.ssfem.N_xi; % Number of stochastic samples
preProc_Variables = BVP.preProc; % Preprocessing variables
E_vector = BVP.ssfem.E_i; % Young's modulus realizations
xs = BVP.ssfem.xs; % Smolyak grid points
ws = BVP.ssfem.ws; % Smolyak weights
M_kappa = BVP.ssfem.M_kappa; % Number of terms in KL expansion
GDof = BVP.preProc.GDof; % Degrees of freedom
P_u = BVP.ssfem.P_u; % Number of polynomial chaos terms
Psi_u = BVP.ssfem.Psi_u; % Polynomial chaos basis functions
activeDofs = BVP.preProc.activeDofs; % Active degrees of freedom
PsiSqNorm = BVP.ssfem.PsiSqNorm; % Norms of the polynomial chaos basis functions
nMC = BVP.preProc.nMC; % Number of Monte Carlo simulations

%% Step 1: Solve FEM for each stochastic sample of Young's modulus
u_pc_realization = zeros(GDof, N_xi); % Initialize displacement realizations matrix
for i = 1:N_xi
    % Solve FEM for each realization of Young’s modulus
    u_pc_realization(:, i) = FEM_1DBar_Tipload_Solver(preProc_Variables, E_vector(i, :));
end

%% Step 2: Polynomial Chaos Expansion (PCE) Basis Evaluation
% Create symbolic and numeric variables for PCE
xi_v = cell(1, M_kappa); % Symbolic variables
xi_n = cell(1, M_kappa); % Numeric variables (samples)
for j = 1:M_kappa
    xi_v{j} = sym(sprintf('xi_%d', j)); % Define symbolic variables
    xi_n{j} = xs(:, j); % Store numerical values from Smolyak grid
end

%% Step 3: Compute PC Expansion Coefficients
u_pc_coef = zeros(GDof, P_u); % Initialize coefficient matrix
psi_j_mat = zeros(N_xi, P_u); % Store evaluated polynomial basis functions

for i = 1:P_u
    psi = Psi_u{i}; % Get the i-th polynomial chaos basis function
    psi_j_mat(:, i) = double(subs(psi, xi_v, xi_n)); % Evaluate basis function
    psi_j = psi_j_mat(:, i); % Assign evaluated basis function
    
    % Compute coefficients using weighted projection
    for j = 1:N_xi
        u_pc_coef(:, i) = u_pc_coef(:, i) + u_pc_realization(:, j) * psi_j(j) * ws(j);
    end
    
    % Normalize the coefficients
    u_pc_coef(:, i) = u_pc_coef(:, i) / PsiSqNorm(i);
end

%% Step 4: Extract Active Degrees of Freedom
u_pc_coef_active = u_pc_coef(activeDofs, :);
mu_u_pc = u_pc_coef(:, 1); % Mean displacement
mu_u_pc_active = mu_u_pc(activeDofs); % Mean displacement for active DOFs

%% Step 5: Compute Covariance of the Displacement Response
C_uu_pc = zeros(GDof); % Initialize covariance matrix
for j = 2:P_u
    C_uu_pc = C_uu_pc + PsiSqNorm(j) * u_pc_coef(:, j) * u_pc_coef(:, j)';
end
% Symmetrize covariance matrix
C_uu_pc = 0.5 * (C_uu_pc + C_uu_pc');
CI_u_pc = sqrt(diag(C_uu_pc)) * 1.96; % 95% confidence interval
C_uu_pc_active = C_uu_pc(2:end, 2:end); % Covariance for active DOFs

%% Step 6: Generate Monte Carlo Samples for Response Displacement
uMC = zeros(GDof, nMC); % Storage for Monte Carlo samples
rng(100, 'twister') % Set random seed for reproducibility
xi = randn(nMC, M_kappa); % Generate standard normal random samples

% Evaluate polynomial basis at Monte Carlo samples
xi_mc = num2cell(xi, 1);
psi_j_ext = zeros(nMC, P_u);
for j = 1:P_u
    psi = Psi_u{j}; % Polynomial basis function
    psi_j_ext(:, j) = double(subs(psi, xi_v, xi_mc)); % Evaluate at Monte Carlo samples
    d_j = u_pc_coef(:, j); % Coefficients from PCE
    for i = 1:GDof
        uMC(i, :) = uMC(i, :) + d_j(i) * psi_j_ext(:, j)'; % Construct displacement field
    end
end

%% Step 7: Compute PDF of Tip Displacement using KDE
uMC_tip = uMC(end, :); % Extract displacement at tip
[pdf_u_tip, x2] = ksdensity(uMC_tip); % Compute probability density function

%% Step 8: Reconstruct Displacement Field
u_xi = u_pc_coef * psi_j_mat.';

%% Assigning results back to BVP structure
BVP.ssfem.mu_u_pc = mu_u_pc;
BVP.ssfem.mu_u_pc_active = mu_u_pc_active;
BVP.ssfem.C_uu_pc = C_uu_pc;
BVP.ssfem.C_uu_pc_active = C_uu_pc_active;
BVP.ssfem.CI_u_pc = CI_u_pc;
BVP.ssfem.u_pc_realization = u_pc_realization;
BVP.ssfem.u_xi = u_xi;
BVP.ssfem.u_pc_coef = u_pc_coef;
BVP.ssfem.u_pc_coef_active = u_pc_coef_active;
BVP.ssfem.pdf_u_tip.pdf = pdf_u_tip;
BVP.ssfem.pdf_u_tip.u = x2;
BVP.statfem.xi = xi;
BVP.ssfem.psi_j_ext = psi_j_ext;
BVP.ssfem.uMC = uMC;

end
