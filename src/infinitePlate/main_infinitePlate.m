%% main_plate_with_hole.m
% --------------------------------------------------------------
% Main script for the 2D Plate with Hole problem using PC-Based-statFEM
%
% Overview:
% This script performs uncertainty quantification and model updating
% for a two-dimensional plate with a hole under load using the 
% Polynomial Chaos-based Statistical Finite Element Method (PC-Based-statFEM).
%
% Steps included:
% 1. Preprocessing: Geometry, material properties, and mesh generation.
% 2. Deterministic FEM Analysis:
%    - Linear Elastic analysis.
%    - Neo-Hookean (nonlinear) analysis.
% 3. KL Expansion: Modeling the spatial variability of Young's modulus.
% 4. Spectral Stochastic FEM (SSFEM) using Polynomial Chaos expansions.
% 5. Model-Reality Mismatch modeling via KL expansion at sensor locations.
% 6. Observational noise modeling using Polynomial Chaos.
% 7. Sensor projection matrix calculation.
% 8. Synthetic observation generation incorporating model mismatch and noise.
% 9. Extended Polynomial Chaos basis construction.
% 10. Hyperparameter estimation via maximum likelihood.
% 11. Bayesian updating of displacement fields with the estimated hyperparameters.
% 12. Visualization of results.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Workspace Cleanup
clearvars();           % Clear all existing variables
clc;                   % Clear the command window
close all;             % Close all open figure windows
addpath('../../lib/'); % Add library path for required functions

% Initialize Boundary Value Problem (BVP) structure
BVP = [];

%% 1. Preprocessing
% Set up geometry, mesh, material properties, and boundary conditions.
BVP = plateWithHole_FEM_preprocess(BVP);

%% 2. Deterministic FEM Analysis
% (a) Linear Elastic FEM Analysis
BVP = plateWithHole_FEM_processLE(BVP);

% (b) Nonlinear FEM Analysis with Neo-Hookean material model
BVP = plateWithHole_FEM_processNH(BVP);

%% 3. KL Expansion of Young's Modulus
% Apply Karhunen-Loeve expansion to model uncertainty in material properties.
BVP = plateWithHole_FEM_processLE_KL_Expansion(BVP);

%% 4. Spectral Stochastic FEM (SSFEM) using Polynomial Chaos expansions
BVP = plateWithHole_SSFEM_process_LE_PC_Expansion(BVP);

%% 5. Model-Reality Mismatch Modeling
% Use KL expansion at sensor locations to account for model-reality discrepancies.
BVP = plateWithHole_statFEM_modelReality_KL_Expansion(BVP, 3);

%% 6. Observational Noise Modeling
% Expand noise using Polynomial Chaos expansions to account for measurement uncertainty.
BVP = plateWithHole_statFEM_Noise_PC_Expansion(BVP);

%% 7. Projection Matrix Calculation
% Compute the matrix H mapping global displacements to sensor measurements.
BVP = plateWithHole_H_matrix(BVP);

%% 8. Synthetic Observation Generation
% Generate synthetic observation data including model-reality mismatch and noise.
BVP = plateWithHole_obs_generate(BVP, 1);

%% 9. Extended Polynomial Chaos Basis Construction
% Extend the PC basis functions to include all uncertainty sources.
BVP = plateWithHole_extended_PC(BVP);

%% 10. Hyperparameter Estimation
% Estimate hyperparameters to align model predictions with observed data.
BVP = plateWithHole_hyperParameter(BVP);

%% 11. Bayesian Updating of Displacement Field
% Update displacement and uncertainty using the estimated hyperparameters.
BVP = plateWithHole_statFEM_PC_update(BVP);

%% 12. Visualization
% Plot the updated displacement field and compare with prior results.
BVP = plateWithHole_FEM_surfplotV2(BVP);

%% Uncomment the following line to calculate and compare the displacement norms:
BVP = plateWithHole_statFEM_norm(BVP);

%% Simulation Completed
disp('Plate with hole statFEM analysis completed!');
