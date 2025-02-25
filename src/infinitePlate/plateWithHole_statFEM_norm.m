function BVP = plateWithHole_statFEM_norm(BVP)
% This function returns the norm of displacement of every
% nodes of the plate from posterior results and true.
% The purpuse of this norm is to show how good the posterior
% displacement field is compared to the true displacement field.
% The idea is originated from the paper:
% "Automated discovery of generalized standard material models with
% EUCLID", Fig. 5
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign variables from BVP
u_truth = BVP.proc.NH.u;
u_post = BVP.statfem.mu_u_y_pc;


%% Calculate the norm
u_truth_x = u_truth(1:2:end);
u_truth_y = u_truth(2:2:end);
u_post_x = u_post{1};
u_post_y = u_post{2};

norm_u_truth = zeros(BVP.preProc.msh.nuNodes,1);
norm_u_post = zeros(BVP.preProc.msh.nuNodes,1);
for inode = 1:BVP.preProc.msh.nuNodes
    norm_u_truth(inode) = norm([u_truth_x(inode) u_truth_y(inode)],2);
    norm_u_post(inode) = norm([u_post_x(inode) u_post_y(inode)],2);
end
norm_u_post_all = norm(norm_u_post)
norm_u_truth_all = norm(norm_u_truth);
%% Save the results
BVP.statfem.norm_u_truth = norm_u_truth;
BVP.statfem.norm_u_post = norm_u_post;

disp('12. Displacement norms are calculated.');
end
