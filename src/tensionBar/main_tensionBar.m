%% main_tension_bar.m
% --------------------------------------------------------------
%  Main script for 1D tension bar under tip load using statFEM
%
% Overview:
% This script performs uncertainty quantification and model updating for a
% one-dimensional tension bar under a tip load using the Statistical Finite
% Element Method (statFEM). It includes:
% - Deterministic FEM analysis
% - Karhunen-Loève (KL) expansion for Young's modulus as weakly-stationary Gaussian random fields
% - Spectral Stochastic FEM (SSFEM) of displacement field via Polynomial Chaos (PC)
% - Model-reality mismatch modeling as non-stationary Gaussian random fields
% - Model-reality mismatch expansion with  Karhunen-Loève (KL)
% - Synthetic observation generation
% - Hyperparameter estimation
% - Bayesian updating of the solution with Gauss-Markov-Kalman filter (GMKF)
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

%% Workspace Cleanup
clearvars();   % Clear all variables
clc;           % Clear command window
close all;     % Close all figure windows
addpath('../../lib/');  % Add library path

% Initialize Boundary Value Problem (BVP) structure
BVP = [];
calcase = 'nonlinear';  % or 'linear'

%% Preprocessing
% Pre-process the model, set geometry, material properties, and boundary conditions
i_ld = 2; % Correlation length index: corrLength_vector = [10, 25, 50, 100]
BVP = tensionBar_FEM_preprocess(BVP, i_ld, calcase);

%% Deterministic FEM Analysis
% Solve the linear elastic problem deterministically
BVP = tensionBar_FEM_process_LE(BVP,calcase);

%% Random Field Modeling (KL Expansion)
% Model spatial variability of Young's Modulus using Karhunen-Loève expansion
BVP = tensionBar_FEM_preprocess_KL_Expansion(BVP);

%% Spectral Stochastic FEM (SSFEM) via Polynomial Chaos
% Incorporate uncertainty in displacement field using Polynomial Chaos expansions
BVP = tensionBar_SSFEM_process_LE_PC_Expansion(BVP);

%% Model-Reality Mismatch Modeling
% Model discrepancy between the computational model and real observations using KL expansion
i_nRed = 2; % Reading level index: 
nRed_vector = [1000, 100, 10, 1];
BVP = tensionBar_statFEM_modelReality_KL_Expansion_converge(BVP, nRed_vector(i_nRed), i_ld, calcase);

%% Observational Noise Modeling
% Expand observational noise using Polynomial Chaos expansions
BVP = tensorBar_statFEM_Noise_PC_Expansion(BVP);

%% Projection Matrix Computation
% Calculate spatial projection matrix (H) to map global displacements to sensor data
BVP = tensionBar_H_matrix(BVP);

%% Synthetic Data Generation
% Generate synthetic observation data incorporating model-reality mismatch and noise
BVP = tensionBar_obs_generate(BVP, i_ld, i_nRed, calcase);

%% Extended Polynomial Chaos Basis
% Extend PC basis functions to capture all sources of uncertainty in a uniform way
BVP = tensionBar_extended_PC(BVP);

%% Hyperparameter Estimation
% Estimate hyperparameters via maximum likelihood to align model predictions with data
BVP = tensionBar_hyperParameter(BVP);

%% Bayesian Update with statFEM
% Perform Bayesian updating of FEM solutions using observation data
BVP = tensionBar_statFEM_PC_update(BVP);

%% Model Evaluation: Root Mean Square Deviation (RMSD)
% Evaluate the updated model accuracy
% BVP = tensionBar_RMSD_nr(BVP, i_ld); % Just for nsen = 9;

% Alternative evaluation method (with sampled data):
% BVP = tensionBar_RMSD_nr_v2(BVP, i_ld); % Just for nsen = 9;
%%
BVP = tensionBar_plot(BVP);
%% Simulation Completed
disp('Tension bar statFEM analysis completed!');
