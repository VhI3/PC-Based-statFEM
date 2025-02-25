%% plateWithHole_extended_PC.m
% -------------------------------------------------------------------------
% Statistical FEM: Extended Polynomial Chaos Expansion (PCE)
%
% Description:
% This function extends the Polynomial Chaos (PC) expansions to incorporate:
%   - Displacement field uncertainties
%   - Model-reality mismatch via the KL expansion
%   - Sensor noise contributions
%
% The extension combines all sources of uncertainty into a unified PC basis.
%
% Inputs:
%   BVP : Boundary Value Problem structure containing:
%         - Existing PC expansions (displacement, sensor noise)
%         - KL expansions (model-reality mismatch)
%         - Observed data and norms
%
% Outputs:
%   BVP : Updated structure containing:
%         - Extended PC basis and coefficients
%         - Updated norms and summed observation projections
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0
% -------------------------------------------------------------------------
function BVP = plateWithHole_extended_PC(BVP)
%% Assign from BVP
% Extracting necessary variables from the BVP structure for extended Polynomial Chaos (PC) expansion
Psi_com = BVP.statfem.Psi_com; % Composite Psi vector from statistical finite element method (STATFEM)
PsiSqNorm = BVP.ssfem.PsiSqNorm; % Square norm of Psi from stochastic spectral finite element method (SSFEM)
ChiSqNorm = BVP.statfem.ChiSqNorm; % Square norm of Chi from sensor noise model
nSen = BVP.preProc.statfem.nuSensors; % Number of sensors
u_pc_coef = BVP.ssfem.u_pc_coef; % Polynomial Chaos coefficients for displacement
u_pc_coef_active = BVP.ssfem.u_pc_coef_active; % Active Polynomial Chaos coefficients for displacement
GDof = BVP.preProc.msh.GDOFs; % Global degrees of freedom
P_u = BVP.ssfem.P_u; % Order of Polynomial Chaos expansion for displacement
e_pc_coef = BVP.statfem.e_pc_coef; % Polynomial Chaos coefficients for noise
M_d = BVP.statfem.M_d; % Number of terms in KL expansion for the sensor
P_d = BVP.statfem.P_d; % Order of Polynomial Chaos expansion for the sensor
sum_yObs = BVP.statfem.sum_yObs; % Sum of observed data
size_activeDOFs = BVP.preProc.BC.size_activeDOFs;
%% Calculation
% Convert Psi composite functions to string and find unique functions for extended PC
Psi_com_str = cellfun(@char, Psi_com, 'UniformOutput', false);
uniq_Psi_com_str = unique(Psi_com_str, 'stable');
Psi_ext = num2cell(str2sym(uniq_Psi_com_str)); % Convert unique strings back to symbolic form
P_ext = length(Psi_ext); % Length of the extended Psi

% Initialize and construct the square norm for the extended Psi
PsiSqNorm_com = []; % Initialize an empty array for composite Psi square norm
PsiSqNorm_ext = [PsiSqNorm_com; PsiSqNorm; ChiSqNorm; ones(nSen, 1)]; % Concatenate square norms for extended Psi

% Extend the u_pc coefficients to match the extended Psi
u_ext_pc_coef = zeros(GDof, P_ext); % Initialize extended u_pc coefficients
u_ext_pc_coef(:, 1:P_u) = u_pc_coef; % Assign existing coefficients


u_ext_pc_cof = {};
u_ext_pc_cof{1} = u_ext_pc_coef(1:2:end,:);
u_ext_pc_cof{2} = u_ext_pc_coef(2:2:end,:);

% Extend the u_pc coefficients for active DOFs
u_ext_pc_coef_active = zeros(size_activeDOFs, P_ext); % Initialize extended active u_pc coefficients
u_ext_pc_coef_active(:, 1:P_u) = u_pc_coef_active; % Assign existing active coefficients

% Extend the e_pc coefficients for noise to match the extended Psi
e_ext_pc_coef = zeros(nSen, P_ext); % Initialize extended e_pc coefficients for noise
e_ext_pc_coef(:, P_u + M_d * P_d + 1:end) = e_pc_coef; % Assign existing noise coefficients

% Initialize and assign the sum of observed data to the extended y_obs PC coefficients
sum_yObs_pc = {};
sum_yObs_pc{1} = zeros(nSen, P_ext); % Initialize sum of observed data in PC form
sum_yObs_pc{2} = zeros(nSen, P_ext); % Initialize sum of observed data in PC form
sum_yObs_pc{1}(:, 1) = sum_yObs{1}; % Assign sum of observed data to the first column
sum_yObs_pc{2}(:, 1) = sum_yObs{2}; % Assign sum of observed data to the first column

%% Assign back to BVP
BVP.statfem.Psi_ext = Psi_ext;
BVP.statfem.PsiSqNorm_ext = PsiSqNorm_ext;
BVP.statfem.P_ext = P_ext;
BVP.statfem.u_ext_pc_coef = u_ext_pc_coef;
BVP.statfem.e_ext_pc_coef = e_ext_pc_coef;
BVP.statfem.u_ext_pc_coef_active = u_ext_pc_coef_active;
BVP.statfem.sum_yObs_pc = sum_yObs_pc;
BVP.statfem.u_ext_pc_cof = u_ext_pc_cof;

disp('9. Extended PC expansion is completed.');
end
