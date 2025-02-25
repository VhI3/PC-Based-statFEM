function BVP = tensionBar_obs_generate(BVP, i_ld, i_nRed, calcase)
% tensionBar_obs_generate - Generates synthetic observation data for the tension bar problem.
%
% This function creates synthetic measurements (observations) of the displacement field
% at sensor locations. The observations include contributions from:
%   - Scaled mean displacement field (from the analytical solution).
%   - Model-reality mismatch (from KL expansion samples).
%   - Measurement noise (sensor error samples).
%
% Inputs:
%   BVP     - Boundary Value Problem structure containing precomputed fields and parameters.
%   i_ld    - Index indicating the correlation length case for saving data.
%   i_nRed  - Index indicating the number of realizations (nRed) for saving data.
%
% Outputs (stored in BVP.statfem):
%   yObs        - Synthetic observation matrix [nSensors x nRed].
%   sum_yObs    - Sum of observations across all realizations [nSensors x 1].
%   mean_yObs   - Mean of the observations [nSensors x 1].
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Extract Variables from BVP Structure

nSen             = BVP.statfem.nSen;          % Number of sensors
nRed             = BVP.statfem.nRed;          % Number of realizations
rho_preDefined   = BVP.preProc.rho_preDefined; % Scaling factor (model bias)
H_full           = BVP.statfem.H_full;        % Projection matrix (global field -> sensor locations)
d_kl_sampled     = BVP.statfem.d_kl_sampled;  % Model-reality mismatch samples (from KL expansion)
e_sampled        = BVP.statfem.e_sampled;     % Measurement noise samples
U_analytic       = BVP.analyticalSolution.U_analytic; % Analytical solution for comparison
U_nonlinear      = BVP.fem.proc.displacement;
%% Generate Synthetic Observations

% Initialize observation matrix: Each column represents a realization.
yObs = zeros(nSen, nRed);

switch calcase
    case 'linear'
        for i = 1:nRed
            % Observation model:
            % yObs = Scaled mean field + Model-reality mismatch + Measurement noise
            yObs(:, i) = rho_preDefined * H_full * U_analytic + d_kl_sampled(:, i) + e_sampled(:, i);
        end
    case 'nonlinear'
        for i = 1:nRed
            % Observation model:
            % yObs = Scaled mean field + Model-reality mismatch + Measurement noise
            yObs(:, i) = rho_preDefined * H_full * U_nonlinear + e_sampled(:, i);
        end
    otherwise
        error('Unknown calculation case type selected.');
end


%% Compute Summary Statistics of Observations

sum_yObs  = sum(yObs, 2);        % Sum of all realizations at each sensor
mean_yObs = sum_yObs / nRed;     % Mean observation at each sensor

%% Save Mean Observations (for specific cases)

if i_nRed == 1
    % Save the mean observation vector based on the correlation length index
    switch i_ld
        case 1
            save mean_yObs_ld10_1000.mat mean_yObs
        case 2
            save mean_yObs_ld25_1000.mat mean_yObs
        case 3
            save mean_yObs_ld50_1000.mat mean_yObs
        case 4
            save mean_yObs_ld100_1000.mat mean_yObs
    end
end

%% Update BVP Structure

BVP.statfem.yObs        = yObs;       % Synthetic observations
BVP.statfem.sum_yObs    = sum_yObs;   % Summation across realizations
BVP.statfem.mean_yObs   = mean_yObs;  % Mean across realizations

end
