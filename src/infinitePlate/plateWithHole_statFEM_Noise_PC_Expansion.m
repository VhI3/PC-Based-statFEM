%% plateWithHole_statFEM_Noise_PC_Expansion.m
% --------------------------------------------------------------
% Statistical FEM: Polynomial Chaos Expansion of Sensor Noise (2D)
%
% Description:
% This function models the sensor measurement noise using a Polynomial Chaos (PC) expansion.
% The sensor noise is assumed to be Gaussian, independent, and identically distributed
% with a known standard deviation. The noise is expanded using Hermite polynomials.
%
% Inputs:
%   BVP : Boundary Value Problem structure with preprocessing data.
%
% Outputs:
%   BVP : Updated structure containing noise PC coefficients, covariance matrix,
%         and sampled sensor errors.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

function BVP = plateWithHole_statFEM_Noise_PC_Expansion(BVP)

%% 1. Extract Variables from BVP
P_e = BVP.preProc.P_e;                           % Order of Hermite Polynomial Chaos (PC) for sensor error
nSen = BVP.preProc.statfem.nuSensors;            % Number of sensors
nRed = BVP.statfem.nRed;                         % Number of reduced samples
sig_e = BVP.preProc.sig_e;                       % Standard deviation of sensor error
Psi_com = BVP.statfem.Psi_com;                   % Composite PC basis functions
nMC = BVP.preProc.nMC;                           % Number of Monte Carlo samples
psi_j_ext = BVP.statfem.psi_j_ext;               % Extended PC basis
sensorNull = BVP.preProc.statfem.sensorNull;     % Indices of inactive sensors (no measurements)

%% 2. Polynomial Chaos Expansion for Sensor Noise
% Generate Hermite PC basis functions for sensor error
[Psi_e] = Hermite_PC_e(nSen, P_e);  

% Error covariance matrix assuming independent and identically distributed noise
e_pc_coef = sig_e * eye(nSen);                  % PC coefficient matrix for sensor noise
C_e_PC = e_pc_coef * e_pc_coef';                % Covariance matrix of sensor noise

% Update composite PC basis with new Psi_e
Psi_com = [Psi_com; Psi_e];

%% 3. Monte Carlo Sampling of Sensor Noise
rng(300, 'twister');                            % Ensure reproducibility

% Generate random samples for sensor noise
zeta_sample = randn(nMC, nSen);                 % Full MC samples for noise
zeta = zeta_sample(1:nRed, :);                  % Reduced set of samples for analysis

% Initialize matrix for sampled noise
e_sampled = sig_e * zeta';                      % Scale samples by standard deviation

% Nullify noise for inactive sensors
e_sampled(sensorNull, :) = 0;

% Extend the PC basis evaluations with the new samples
psi_j_ext = [psi_j_ext, zeta_sample];

%% 4. Assign Results to BVP Structure
BVP.statfem.Psi_e = Psi_e;                      % Hermite PC basis for noise
BVP.statfem.e_pc_coef = e_pc_coef;              % Noise PC coefficients
BVP.statfem.C_e_PC = C_e_PC;                    % Noise covariance matrix
BVP.statfem.e_sampled = e_sampled;              % Sampled sensor noise
BVP.statfem.Psi_com = Psi_com;                  % Updated composite PC basis
BVP.statfem.psi_j_ext = psi_j_ext;              % Updated extended PC evaluations

disp('7. Sensor noise PC expansion completed.');

end

