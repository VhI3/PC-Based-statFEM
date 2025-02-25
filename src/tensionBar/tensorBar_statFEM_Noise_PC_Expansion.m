function BVP = tensorBar_statFEM_Noise_PC_Expansion(BVP)
% This function performs the Polynomial Chaos (PC) expansion of the sensor error 
% (measurement noise) in the statFEM framework. The sensor error is modeled as a 
% Gaussian random field with a given standard deviation and expanded using Hermite 
% polynomials to incorporate uncertainty in the measurement process.
%
% Inputs:
%   BVP - Boundary Value Problem structure containing necessary preprocessing 
%         information and previously computed variables.
%
% Outputs:
%   Updated BVP structure with:
%       - Polynomial Chaos basis functions for sensor error (Psi_e)
%       - PC coefficients for the error (e_pc_coef)
%       - Error covariance matrix (C_e_PC)
%       - Sampled sensor errors (e_sampled)
%       - Updated composite Psi matrix (Psi_com)
%       - Extended psi_j_ext with sampled variables
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Extract Required Variables from BVP
P_e = BVP.preProc.P_e;          % Order of Hermite Polynomial Chaos (PC) for sensor error
nSen = BVP.statfem.nSen;        % Number of sensors
nRed = BVP.statfem.nRed;        % Number of reduced random samples
sig_e = BVP.preProc.sig_e;      % Standard deviation of the sensor error
Psi_com = BVP.statfem.Psi_com;  % Composite Psi matrix accumulated from previous expansions
nMC = BVP.preProc.nMC;          % Number of Monte Carlo simulations
psi_j_ext = BVP.statfem.psi_j_ext; % Extended basis functions from previous steps

%% Generate Hermite Polynomial Chaos Basis for Sensor Error
% Psi_e: Hermite polynomial basis functions used to expand the sensor error.
[Psi_e] = Hermite_PC_e(nSen, P_e);  

%% Construct the Error Covariance Matrix
% Assuming independent sensor errors with equal variance:
e_pc_coef = sig_e * eye(nSen);     % Coefficient matrix: scaled identity (diagonal covariance)
C_e_PC = e_pc_coef * e_pc_coef';   % Covariance matrix: diagonal with sig_e^2 on the diagonal

%% Update Composite Psi Matrix
% Append the new Psi_e matrix (sensor error basis) to the overall composite matrix.
Psi_com = [Psi_com; Psi_e];       

%% Generate Random Samples for Sensor Error
rng(300, 'twister');               % Set random seed for reproducibility

% zeta: random samples for reduced order modeling
% zeta_sample: random samples for Monte Carlo simulations
zeta = randn(nRed, nSen);          
zeta_sample = randn(nMC, nSen);    

%% Compute Sampled Sensor Errors
% e_sampled: matrix where each column represents a sampled sensor error realization.
e_sampled = zeros(nSen, nRed);     
for s = 1:nSen
  e_sampled(s, :) = sig_e * zeta(:, s);  % Scale standard normal samples by sig_e
end

% Uncomment to calculate the sample covariance matrix from the generated samples:
% C_e_sampled = cov(e_sampled);

%% Extend Polynomial Chaos Evaluations
% Extend the existing Psi evaluations with the newly sampled sensor error variables.
psi_j_ext = [psi_j_ext, zeta_sample];    

%% Assign Updated Variables Back to BVP Structure
BVP.statfem.Psi_e = Psi_e;           % Hermite PC basis functions for sensor error
BVP.statfem.e_pc_coef = e_pc_coef;   % Error PC coefficients
BVP.statfem.C_e_PC = C_e_PC;         % Covariance matrix of the sensor error
BVP.statfem.e_sampled = e_sampled;   % Sampled sensor errors
BVP.statfem.Psi_com = Psi_com;       % Updated composite Psi matrix
BVP.statfem.zeta = zeta;             % Reduced random variables for sensor noise
BVP.statfem.psi_j_ext = psi_j_ext;   % Updated extended basis evaluations

end
