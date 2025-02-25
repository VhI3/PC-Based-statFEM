%% plateWithHole_statFEM_modelReality_KL_Expansion.m
% --------------------------------------------------------------
% Statistical FEM: Model Reality Mismatch via KL Expansion (2D)
%
% Description:
% This function models the discrepancy between the computational model and the
% physical reality using a Karhunen-Loeve (KL) expansion at sensor locations.
% It uses a Fredholm-Galerkin method to compute the eigenfunctions and eigenvalues
% of the covariance kernel and represents the stochastic field with polynomial chaos.
%
% Inputs:
%   BVP       : Boundary Value Problem structure with preprocessing data.
%   obserCase : Observation case controlling sensor configuration and sample sizes.
%
% Outputs:
%   BVP : Updated structure with KL expansion parameters, sensor eigenvalues,
%         eigenvectors, and extended PC basis functions.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

function BVP = plateWithHole_statFEM_modelReality_KL_Expansion(BVP,obserCase)
%% Assign from BVP
preProc_Variables                    = BVP.preProc;
L                                    = preProc_Variables.geometry.L;
H                                    = preProc_Variables.geometry.H;
Psi_com = BVP.ssfem.Psi_com; % Composite Psi matrix from STATFEM
nMC = BVP.preProc.nMC;
psi_j_ext = BVP.statfem.psi_j_ext;
P_d = 1;
%
%
%
%% Karhunen-Loeve Expansion and Polynomial Chaos Calculations
% Perform Karhunen-Loeve expansion using the Fredholm-Galerkin method to obtain
% numerical eigenvalues and eigenvectors for the expansion.
% Setting parameters based on observation case
% This switch case sets the reduction factor, number of sensors, dot size, and alpha value based on the observation case
switch obserCase
  case 1
    nRed = 1;  dotSize = 30; allp = 0.1;
  case 2
    nRed = 11;  dotSize = 15; allp = 0.1;
  case 3
    nRed = 100;  dotSize = 5; allp = 0.1;
  case 4
    nRed = 200;  dotSize = 1; allp = 0.1;
end

Lx_d            = L/4;
Ly_d            = H/4;


cov_func_d = @(x1, x2) (1 + sqrt(5)*abs(x1(1) - x2(1))/Lx_d + 5/3*(abs(x1(1) - x2(1))^2)/(Lx_d ^ 2))*...
    (1 + sqrt(5)*abs(x1(2) - x2(2))/Ly_d + 5/3*(abs(x1(2) - x2(2))^2)/(Ly_d ^ 2))*...
    exp(-sqrt(5)*abs(x1(1) - x2(1))/Lx_d - sqrt(5)*abs(x1(2) - x2(2))/Ly_d);


preProc_Variables_d.msh.nElm             = BVP.preProc.statfem.nuSenElm;
preProc_Variables_d.msh.nodeCoordinates      = BVP.preProc.statfem.sensorCoordinates;
preProc_Variables_d.msh.elementNodes              = BVP.preProc.statfem.sensorNodes;
preProc_Variables_d.proc.nuNodes         = BVP.preProc.statfem.nuSensors;
preProc_Variables_d.proc.order2D         = preProc_Variables.proc.order2D;

tic
[sensor_eigVal, ~, sensor_eigenVec, sensor_MLterms] = ...
  KL_fredholm_Galerkin_local_2D_optimized(preProc_Variables_d, cov_func_d);
toc

M_d = sensor_MLterms;

% Initialize and compute the square root of eigenvalues times eigenvectors
sensor_eigMat = zeros(size(sensor_eigenVec));
for i = 1:M_d
  sensor_eigMat(:, i) = sensor_eigenVec(:, i) * sqrt(sensor_eigVal(i));
end


[Psi_d, ~] = Hermite_PC_d(M_d, P_d);


rng(200,'twister');
chi = randn(nRed,M_d);
chi_sample = randn(nMC,M_d);
psi_d = cell(P_d, M_d); % Initialize cell array for Hermite PC evaluations

% Evaluate the Hermite PC basis functions at the sampled points
for k = 1:M_d
  chi_n = chi(:, k); % Sample points for the k-th term
  chi_v = sym(sprintf('chi_%d', k)); % Symbolic variable for the k-th term
  chi_sample_n = chi_sample(:, k);
  for l = 1:P_d
    psi = Psi_d{l, k}; % Hermite PC basis function
    psi_d{l, k} = double(subs(psi, chi_v, chi_n)); % Evaluate basis function at sample points
    psi_j_ext = [psi_j_ext double(subs(psi, chi_v, chi_sample_n))];
  end
end



% Compute the square norm matrix for the Chi variables
ChiSqNorm_matrix = zeros(P_d * M_d);
counter = 1;
for iM_d = 1:M_d
  for iP_d = 1:P_d
    ChiSqNorm_matrix(counter, counter) = factorial(iP_d); % Factorial for normalization
    counter = counter + 1;
  end
end
ChiSqNorm = diag(ChiSqNorm_matrix);

% Append the reshaped Psi_d to the composite Psi matrix
Psi_com = [Psi_com; reshape(Psi_d, M_d * P_d, 1)];


%% Assign back to BVP
BVP.statfem.nRed = nRed;
BVP.statfem.M_d = M_d;
BVP.statfem.P_d = P_d;
BVP.statfem.sensor_eigVal = sensor_eigVal;
BVP.statfem.sensor_eigenVec = sensor_eigenVec;
BVP.statfem.sensor_eigMat = sensor_eigMat;
BVP.statfem.Psi_d = Psi_d;
BVP.statfem.psi_j_ext = psi_j_ext;
BVP.statfem.Psi_com = Psi_com;
BVP.statfem.ChiSqNorm_matrix = ChiSqNorm_matrix;
BVP.statfem.ChiSqNorm = ChiSqNorm;
BVP.plotSettings.dotSize = dotSize;
BVP.plotSettings.allp = allp;

disp('6. KL expansion of model discrepancy is completed.');
end

