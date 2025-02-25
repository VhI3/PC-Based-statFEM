function BVP = tensionBar_RMSD_nr(BVP, i_ld)
% This function calculates the Root Mean Square Deviation (RMSD) between the posterior mean
% displacement at the sensor locations (mu_z_pc) and the reference mean observation (mean_yObs).
% RMSD is used to assess the accuracy of the updated model with respect to synthetic observation data.

%% Assign from BVP
% Load the pre-saved mean observation data based on the correlation length index (i_ld)
if i_ld == 1
    load mean_yObs_ld10_1000.mat   % For correlation length ld = 10
elseif i_ld == 2
    load mean_yObs_ld25_1000.mat   % For correlation length ld = 25
elseif i_ld == 3
    load mean_yObs_ld50_1000.mat   % For correlation length ld = 50
elseif i_ld == 4
    load mean_yObs_ld100_1000.mat  % For correlation length ld = 100
end

% Extract posterior mean displacement at sensor locations from BVP
mu_z_pc = BVP.statfem.mu_z_pc;  

%% Calculate RMSD
% The RMSD is computed as the Euclidean norm of the difference between the reference mean observation 
% (mean_yObs) and the posterior mean prediction (mu_z_pc).

% RMSD calculation (difference norm between posterior prediction and reference data)
RMSE_ = norm(mean_yObs - mu_z_pc);  

%% Assign back to BVP
BVP.RMSE_ = RMSE_;  % Store calculated RMSD in BVP structure
end
