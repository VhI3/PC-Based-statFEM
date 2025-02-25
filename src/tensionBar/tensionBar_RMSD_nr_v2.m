function BVP = tensionBar_RMSD_nr_v2(BVP, i_ld)
% This function calculates the Root Mean Square Deviation (RMSD) between the Monte Carlo (MC)
% simulated posterior displacement field and the reference mean observation data.
% RMSD is used to quantify how well the updated model predictions match the synthetic observations.

%% Assign from BVP
% Load pre-saved reference mean observation data based on the correlation length index (i_ld)
if i_ld == 1
    load mean_yObs_ld10_1000.mat   % For correlation length ld = 10
elseif i_ld == 2
    load mean_yObs_ld25_1000.mat   % For correlation length ld = 25
elseif i_ld == 3
    load mean_yObs_ld50_1000.mat   % For correlation length ld = 50
elseif i_ld == 4
    load mean_yObs_ld100_1000.mat  % For correlation length ld = 100
end

% Extract necessary variables from the BVP structure
nMC = BVP.preProc.nMC;            % Number of Monte Carlo samples
u_y_MC = BVP.statfem.u_y_MC;      % Posterior displacement samples from Monte Carlo simulation (size: GDof x nMC)
H_full = BVP.statfem.H_full;      % Projection matrix to map global displacements to sensor locations

%% Calculate RMSD
% Initialize sum of squared differences
RMSE_sum = 0;  

% Loop through each Monte Carlo sample to calculate the squared difference
for inMC = 1:nMC
    % Project the displacement sample to sensor locations
    projected_displacement = H_full * u_y_MC(:, inMC);  
    
    % Calculate the squared difference between the projected displacement and the mean observed data
    squared_difference = (projected_displacement - mean_yObs)' * (projected_displacement - mean_yObs);  
    
    % Accumulate the squared differences
    RMSE_sum = RMSE_sum + squared_difference;  
end

% Calculate the mean of the accumulated squared differences
mean_RMSE_sum = RMSE_sum / nMC;  

% Compute the final RMSD by taking the square root
RMSE_ = sqrt(mean_RMSE_sum);  

%% Assign back to BVP
% Store the calculated RMSD in the BVP structure
BVP.RMSE_ = RMSE_;  
end
