%% plateWithHole_obs_generate.m
% --------------------------------------------------------------
% Statistical FEM: Synthetic Observation Generation
%
% Description:
% This function generates synthetic observation data based on the displacement
% field obtained from deterministic FEM solutions. It incorporates both:
%   - Model-reality mismatch (from the chosen option: Linear Elastic (LE) or Neo-Hookean (NH))
%   - Observational noise (sensor error)
%
% Inputs:
%   BVP    : Boundary Value Problem structure containing FEM solutions, sensor info, and noise samples.
%   option : Integer flag specifying the FEM solution type:
%            1 - Linear Elastic (LE) solution with inhomogeneous Young's modulus
%            2 - Neo-Hookean (NH) solution
%
% Outputs:
%   BVP : Updated structure with:
%         - Synthetic observations (X and Y directions)
%         - Observations mapped onto surface mesh
%         - Summed observations over samples
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

function BVP = plateWithHole_obs_generate(BVP,option)
%% Assign from BVP
if option == 1
    ux               = BVP.proc.LE.ux_nonHom;
    uy               = BVP.proc.LE.uy_nonHom;
    u_nonHom         = BVP.proc.LE.u_nonHom;
elseif option ==2
    ux               = BVP.proc.NH.ux;
    uy               = BVP.proc.NH.uy;
    u_nonHom         = BVP.proc.NH.u;
else
    error('Option is not valid');
end
nSam             = BVP.statfem.nRed;
e_sample         = BVP.statfem.e_sampled;
nSen             = BVP.preProc.statfem.nuSensors;
Hx_full           = BVP.obs.Hx_full;
Hy_full          = BVP.obs.Hy_full;
H_full            = BVP.obs.H_full;
sensorNodes       = BVP.preProc.statfem.sensorNodes;


%% Calculation
Y_exp_x = zeros(nSen,nSam);
Y_exp_y = zeros(nSen,nSam);

Y_exp0 = zeros(2*nSen,nSam);
e_sample0  = reshape(repmat(e_sample', 2, 1), size(e_sample, 2), [])';
for i = 1:nSam
    Y_exp_x(:,i) = Hx_full*ux + e_sample(:,i);
    Y_exp_y(:,i) = Hy_full*uy + e_sample(:,i);
    Y_exp0(:,i) =  H_full*u_nonHom + e_sample0(:,i);
    
end
Y_exp = zeros(2*nSen,nSam);
Y_exp(1:2:end,:) = Y_exp_x;
Y_exp(2:2:end,:) = Y_exp_y;

Y_exp_x_Surf = makeSurf(sensorNodes,Y_exp_x(:,1));
Y_exp_y_Surf = makeSurf(sensorNodes,Y_exp_y(:,1));

% Compute the sum and mean of the generated observations across all realizations
sum_yObs = {};
sum_yObs{1} = sum(Y_exp_x, 2);
sum_yObs{2} = sum(Y_exp_y, 2);



%% Assign back to BVP
BVP.statfem.Y_exp_x = Y_exp_x;
BVP.statfem.Y_exp_y = Y_exp_y;
BVP.statfem.Y_exp   = Y_exp;
BVP.statfem.sum_yObs = sum_yObs; % Sum of observations
BVP.statfem.Y_exp_x_Surf = Y_exp_x_Surf;
BVP.statfem.Y_exp_y_Surf = Y_exp_y_Surf;


disp('8. Synthetic observation generation is completed.');
end
