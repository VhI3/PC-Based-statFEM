function BVP = tensionBar_plot(BVP)
% tensionBar_plot: Visualizes the displacement fields and observation data
% for the 1D tension bar problem.
%
% Inputs:
%   - BVP: Boundary value problem structure containing solution data.
%   - plotS: Plot settings (unused in this function but reserved for future).
%   - cc: Color code (unused but reserved for customization).
%
% Outputs:
%   - BVP: The input structure is returned unchanged (for consistency).
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)


%% Extract variables from BVP
nodeCoordinates = BVP.preProc.nodeCoordinates;     % Node coordinates along the bar
nRed = BVP.statfem.nRed;                            % Number of realizations for observation data
mu_u_pc = BVP.ssfem.mu_u_pc;                        % Prior mean displacement field
CI_u_pc = BVP.ssfem.CI_u_pc;                        % 95% confidence interval for prior displacement
yObs = BVP.statfem.yObs;                            % Observed displacement data at sensors
mu_u_y_pc = BVP.statfem.mu_u_y_pc;                  % Posterior mean displacement field after updating
CI_u_y_pc = BVP.statfem.CI_u_y_pc;                  % 95% confidence interval for posterior displacement
L = BVP.preProc.geometry.L;                         % Length of the tension bar
close all;
%% Plotting Prior Displacement Field with Confidence Interval
% Plot the prior mean displacement field (before Bayesian update)
plot_u_FEM_mean = plot(nodeCoordinates, mu_u_pc, 'b-', 'LineWidth', 0.5);
set(plot_u_FEM_mean, 'DisplayName', '$\mu_{u_f}$'); % Legend entry for prior mean
hold on;

% Fill area for the 95% confidence interval of the prior displacement
opts = {'EdgeColor', [0 0 1], 'FaceColor', [0 0 1], 'alpha', 0.5};
valid_range = nodeCoordinates > -1 & nodeCoordinates < nodeCoordinates(end) + 1; % Range for plotting
plot_ci_c_u = fill_between(nodeCoordinates, (mu_u_pc - CI_u_pc), (mu_u_pc + CI_u_pc), valid_range, opts{:});
set(plot_ci_c_u, 'DisplayName', '$95$\% CI of $u_f$'); % Legend entry for prior CI
grid on;

%% Plotting Observed Data at Sensor Locations
% Scatter plot of observed data points
for i = 1:nRed
    h = scatter(BVP.statfem.senCoordinates, yObs(:, i), 'k.', 'SizeData', BVP.plotSettings.dotSize);
    if i ~= nRed
        h.Annotation.LegendInformation.IconDisplayStyle = 'off'; % Show legend only for the last point
    else
        set(h, 'DisplayName', 'obs. data'); % Legend entry for observation data
    end
end

%% Plotting Posterior Displacement Field with Confidence Interval
% Plot the posterior mean displacement field (after Bayesian update)
plot_u_FEM_mean_kalman = plot(nodeCoordinates, mu_u_y_pc, 'r-', 'LineWidth', 0.5);
set(plot_u_FEM_mean_kalman, 'DisplayName', '$\mu_{u_a}$'); % Legend entry for posterior mean

% Fill area for the 95% confidence interval of the posterior displacement
opts = {'EdgeColor', [1 0 0], 'FaceColor', [1 0 0], 'alpha', 0.3};
plot_cu_post_kalman = fill_between(nodeCoordinates, ...
    (mu_u_y_pc - CI_u_y_pc), (mu_u_y_pc + CI_u_y_pc), valid_range, opts{:});
set(plot_cu_post_kalman, 'DisplayName', '$95$\% CI of $u_a$'); % Legend entry for posterior CI

%
opts={'EdgeColor', [0 0 0], 'FaceColor', [0 0 0], 'alpha',0.3};
where = BVP.statfem.senCoordinates>-1 & BVP.statfem.senCoordinates<(BVP.statfem.senCoordinates(end)+1);
plot_cu_z_kalman = fill_between(BVP.statfem.senCoordinates, ...
    (BVP.statfem.mu_z_pc - BVP.statfem.CI_z_pc), (BVP.statfem.mu_z_pc + BVP.statfem.CI_z_pc), where, opts{:});
set(plot_cu_z_kalman,'DisplayName','$95$\% CI of $z$')

hl = legend('show', 'FontSize', 8, 'location', 'northwest');
set(hl, 'Interpreter', 'latex');
%% Axis Labels and Plot Formatting
xlabel('$X$', 'interpreter', 'latex');         % X-axis label
ylabel('$u_f$, $u_a$', 'interpreter', 'latex'); % Y-axis label
max_displacement = max(mu_u_pc);                              % Maximum value for setting plot limits
ylim([-0.1 * max_displacement 40]);                           % Y-axis limits
xlim([-5 L + 5]);                                             % X-axis limits
axis square;                                                  % Make plot square-shaped

end
