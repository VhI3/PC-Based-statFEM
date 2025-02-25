function BVP = tensionBar_extended_PC(BVP)
% tensionBar_extended_PC - Constructs the extended Polynomial Chaos (PC) basis and coefficients.
%
% This function extends the PC basis functions and associated coefficients to include:
%   - Displacement field expansions.
%   - Measurement noise expansions.
%   - Model-reality mismatch expansions.
%
% The extended basis allows consistent handling of all uncertainty sources in a unified PC framework.
%
% Inputs:
%   BVP - Boundary Value Problem structure containing precomputed fields and coefficients.
%
% Outputs (updated in BVP.statfem):
%   Psi_ext              - Extended set of PC basis functions.
%   PsiSqNorm_ext        - Extended square norms of the PC basis functions.
%   P_ext                - Total number of extended PC terms.
%   u_ext_pc_coef        - Extended PC coefficients for the displacement field.
%   u_ext_pc_coef_active - Extended PC coefficients for active degrees of freedom (DOFs).
%   e_ext_pc_coef        - Extended PC coefficients for the measurement noise.
%   sum_yObs_pc          - Extended PC representation of the sum of observed data.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Extract Variables from BVP Structure

Psi_com              = BVP.statfem.Psi_com;           % Composite Psi functions from statFEM (symbolic form)
PsiSqNorm            = BVP.ssfem.PsiSqNorm;           % Square norm of Psi (for displacement PC terms)
ChiSqNorm            = BVP.statfem.ChiSqNorm;         % Square norm of Chi (for noise PC terms)
nSen                 = BVP.statfem.nSen;              % Number of sensors
u_pc_coef            = BVP.ssfem.u_pc_coef;           % PC coefficients for displacement field
u_pc_coef_active     = BVP.ssfem.u_pc_coef_active;    % Active displacement PC coefficients (excluding fixed nodes)
GDof                 = BVP.preProc.GDof;              % Global degrees of freedom
P_u                  = BVP.ssfem.P_u;                 % Number of displacement PC terms
e_pc_coef            = BVP.statfem.e_pc_coef;         % PC coefficients for sensor noise
M_d                  = BVP.statfem.M_d;               % Number of KL terms for model-reality mismatch
P_d                  = BVP.preProc.P_d;               % Number of PC terms for model-reality mismatch
sum_yObs             = BVP.statfem.sum_yObs;          % Sum of observed data across all realizations

%% Process and Extend Polynomial Chaos Basis Functions

% Convert Psi composite functions into string format for comparison
Psi_com_str = cellfun(@char, Psi_com, 'UniformOutput', false);

% Identify unique basis functions (avoiding duplicates)
uniq_Psi_com_str = unique(Psi_com_str, 'stable');

% Convert the unique string representations back to symbolic form
Psi_ext = num2cell(str2sym(uniq_Psi_com_str));

% Total number of extended PC terms
P_ext = length(Psi_ext);

%% Extend the Square Norms for the Combined Basis

% Initialize and combine square norms:
%   - PsiSqNorm: from displacement PC.
%   - ChiSqNorm: from noise PC.
%   - ones(nSen, 1): for direct sensor measurements without random components.
PsiSqNorm_ext = [PsiSqNorm; ChiSqNorm; ones(nSen, 1)];

%% Extend the Displacement PC Coefficients

% Initialize extended PC coefficients for displacement
u_ext_pc_coef = zeros(GDof, P_ext);

% Fill in the existing displacement PC coefficients in the extended matrix
u_ext_pc_coef(:, 1:P_u) = u_pc_coef;

%% Extend the Active Displacement PC Coefficients

% Initialize extended PC coefficients for active DOFs (excluding fixed nodes)
u_ext_pc_coef_active = zeros(GDof - 1, P_ext);

% Populate with the previously calculated active displacement coefficients
u_ext_pc_coef_active(:, 1:P_u) = u_pc_coef_active;

%% Extend the Noise PC Coefficients

% Initialize extended noise PC coefficient matrix
e_ext_pc_coef = zeros(nSen, P_ext);

% Place existing noise coefficients at their appropriate location in the extended basis
% The noise terms start after displacement and model-reality mismatch terms.
e_ext_pc_coef(:, P_u + M_d * P_d + 1:end) = e_pc_coef;

%% Extend the Observed Data into PC Form

% Initialize the extended PC coefficients for the sum of observed data
sum_yObs_pc = zeros(nSen, P_ext);

% Assign the sum of observations to the first PC term (mean component)
sum_yObs_pc(:, 1) = sum_yObs;

%% Update BVP Structure with Extended PC Data

BVP.statfem.Psi_ext                = Psi_ext;              % Extended PC basis
BVP.statfem.PsiSqNorm_ext          = PsiSqNorm_ext;        % Extended square norms
BVP.statfem.P_ext                  = P_ext;                % Total number of extended terms
BVP.statfem.u_ext_pc_coef          = u_ext_pc_coef;        % Extended displacement PC coefficients
BVP.statfem.u_ext_pc_coef_active   = u_ext_pc_coef_active; % Extended active displacement PC coefficients
BVP.statfem.e_ext_pc_coef          = e_ext_pc_coef;        % Extended noise PC coefficients
BVP.statfem.sum_yObs_pc            = sum_yObs_pc;          % Extended PC representation of observation sums

end
