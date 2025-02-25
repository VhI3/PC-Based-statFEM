function BVP = plateWithHole_FEM_surfplotV2(BVP)
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
% Plot patches for the six layers in the mesh
plot_element0_mesh = patch(BVP.preProc.msh.xSurf, BVP.preProc.msh.ySurf, 'black', 'FaceAlpha', .01, 'EdgeColor', 'blue', 'LineWidth', .2, 'LineStyle', '-');
set(plot_element0_mesh, 'DisplayName', 'Elements with $E = 200\, $ GPa');
plot_element5_mesh = patch(BVP.preProc.msh.xSurf_5, BVP.preProc.msh.ySurf_5, 'black', 'FaceAlpha', .1, 'EdgeColor', 'blue', 'LineWidth', .2, 'LineStyle', '-');
set(plot_element5_mesh, 'DisplayName', 'Elements with $E = 100\, $ GPa');
plot_element4_mesh = patch(BVP.preProc.msh.xSurf_4, BVP.preProc.msh.ySurf_4, 'black', 'FaceAlpha', .2, 'EdgeColor', 'blue', 'LineWidth', .2, 'LineStyle', '-');
set(plot_element4_mesh, 'DisplayName', 'Elements with $E = 80\, $ GPa');
plot_element3_mesh = patch(BVP.preProc.msh.xSurf_3, BVP.preProc.msh.ySurf_3, 'black', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'LineWidth', .2, 'LineStyle', '-');
set(plot_element3_mesh, 'DisplayName', 'Elements with $E = 50\, $ GPa');
plot_element2_mesh = patch(BVP.preProc.msh.xSurf_2, BVP.preProc.msh.ySurf_2, 'black', 'FaceAlpha', .4, 'EdgeColor', 'blue', 'LineWidth', .2, 'LineStyle', '-');
set(plot_element2_mesh, 'DisplayName', 'Elements with $E = 20\, $ GPa');
plot_element1_mesh = patch(BVP.preProc.msh.xSurf_1, BVP.preProc.msh.ySurf_1, 'black', 'FaceAlpha', .5, 'EdgeColor', 'blue', 'LineWidth', .2, 'LineStyle', '-');
set(plot_element1_mesh, 'DisplayName', 'Elements with $E = 10\, $ GPa');

box on;
axis equal;
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot_fem_mesh = ...
    patch(BVP.preProc.msh.xSurf, BVP.preProc.msh.ySurf, ...
    'blue','EdgeColor','blue','FaceAlpha',.9,...
    'FaceColor','none','LineWidth',.2,'LineStyle','-.'); %#ok<NASGU>
hold on
plot_homogeneous_deformed = patch(BVP.preProc.msh.xSurf + BVP.proc.LE.ux_Surf, BVP.preProc.msh.ySurf + BVP.proc.LE.uy_Surf, ...
    'blue', 'EdgeColor', 'blue', 'FaceColor', 'none', 'LineWidth', 0.6, 'LineStyle', '-'); %#ok<NASGU>

legend_entries = {'Undeformed Shape', ...
    'Homogeneous Def. Shape with LE'};

% Plot dummy rectangles for legend
h = zeros(length(legend_entries), 1);
h(1) = plot([0 0], [0 0], 'LineWidth', 0.2, 'LineStyle', '-.', 'Color', 'blue');
h(2) = plot([0 0], [0 0], 'LineWidth', 0.6, 'LineStyle', '-' , 'Color', 'blue');


% Show legend
legend(h, legend_entries, 'location', 'northwest', 'Interpreter', 'latex');
box on;
axis equal;
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot_homogeneous_deformed = patch(BVP.preProc.msh.xSurf + BVP.proc.LE.ux_Surf, BVP.preProc.msh.ySurf + BVP.proc.LE.uy_Surf, ...
    'blue', 'EdgeColor', 'blue', 'FaceColor', 'none', 'LineWidth', 0.2, 'LineStyle', '-'); %#ok<NASGU>
hold on
plot_inhomogeneous_deformed = ...
    patch(BVP.preProc.msh.xSurf + BVP.proc.NH.ux_Surf, BVP.preProc.msh.ySurf + BVP.proc.NH.uy_Surf, ...
    'green', 'EdgeColor', [0 0.5 0], 'FaceColor', 'none', 'LineWidth', 0.6, 'LineStyle', '-'); %#ok<NASGU>

legend_entries = {'Homogeneous Def. Shape with LE', ...
    'Inhomogeneous Def. Shape with NH'};

% Plot dummy rectangles for legend
h = zeros(length(legend_entries), 1);
h(1) = plot([0 0], [0 0], 'LineWidth', 0.2, 'LineStyle', '-', 'Color', 'blue');
h(2) = plot([0 0], [0 0], 'LineWidth', 0.6, 'LineStyle', '-', 'Color', [0 0.5 0]);
% Show legend
legend(h, legend_entries, 'location', 'northwest', 'Interpreter', 'latex');
box on;
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
plot_homogeneous_dispX_contour = ...
    patch(BVP.preProc.msh.xSurf,BVP.preProc.msh.ySurf,...
    BVP.proc.LE.ux_Surf,'FaceColor','interp'); %#ok<NASGU>
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
plot_inhomogeneous_dispX_contour = ...
    patch(BVP.preProc.msh.xSurf,BVP.preProc.msh.ySurf,...
    BVP.proc.NH.ux_Surf,'FaceColor','interp'); %#ok<NASGU>
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(51)
plot_observation_dispX_contour = ...
    patch(BVP.preProc.msh.xSurf,BVP.preProc.msh.ySurf,...
    BVP.proc.NH.ux_Surf,'FaceColor','interp'); %#ok<NASGU>
hold on
scatter(BVP.preProc.statfem.sensorCoordinates(:,1), BVP.preProc.statfem.sensorCoordinates(:,2), ...
    'SizeData',5,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',1);
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
% axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
plot_mean_dispX_contour = ...
    patch(BVP.preProc.msh.xSurf + BVP.proc.LE.ux_Surf,BVP.preProc.msh.ySurf + BVP.proc.LE.uy_Surf, ...
    BVP.ssfem.mu_ux_pc_Surf,'FaceColor','interp'); %#ok<NASGU>
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
% % axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)
plot_std_dispX_contour = ...
    patch(BVP.preProc.msh.xSurf + BVP.proc.LE.ux_Surf,BVP.preProc.msh.ySurf + BVP.proc.LE.uy_Surf, ...
    BVP.ssfem.std_ux_pc_Surf,'FaceColor','interp'); %#ok<NASGU>
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
% % axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8)
plot_prior_deformed = patch(BVP.preProc.msh.xSurf + BVP.proc.LE.ux_Surf, BVP.preProc.msh.ySurf + BVP.proc.LE.uy_Surf, ...
    'blue', 'EdgeColor', 'blue', 'FaceColor', 'none', 'LineWidth', 0.2, 'LineStyle', '-'); %#ok<NASGU>
hold on
plot_deformed_Y = ...
    patch(BVP.preProc.msh.xSurf + BVP.statfem.mu_ux_y_pc_Surf, ...
    BVP.preProc.msh.ySurf + BVP.statfem.mu_uy_y_pc_Surf,...
    'red','EdgeColor', 'red', 'FaceColor','none','LineWidth',0.6,'LineStyle','-'); %#ok<NASGU>

legend_entries = {'Prior Deformed Shape', ...
    'Posterior Deformed Shape'};

% Plot dummy rectangles for legend
h = zeros(length(legend_entries), 1);
h(1) = plot([0 0], [0 0], 'LineWidth', 0.2, 'LineStyle', '-', 'Color', 'blue');
h(2) = plot([0 0], [0 0], 'LineWidth', 0.6, 'LineStyle', '-' , 'Color', 'red');


% Show legend
legend(h, legend_entries, 'location', 'northwest', 'Interpreter', 'latex');

% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal;
box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
plot_inhomogeneous_deformed = ...
    patch(BVP.preProc.msh.xSurf + BVP.proc.NH.ux_Surf, BVP.preProc.msh.ySurf + BVP.proc.NH.uy_Surf, ...
    'green', 'EdgeColor', [0 0.5 0], 'FaceColor', 'none', 'LineWidth', 0.3, 'LineStyle', '-'); %#ok<NASGU>
hold on
plot_deformed_Y = ...
    patch(BVP.preProc.msh.xSurf + BVP.statfem.mu_ux_y_pc_Surf, ...
    BVP.preProc.msh.ySurf + BVP.statfem.mu_uy_y_pc_Surf,...
    'red','EdgeColor', 'red', 'FaceColor','none','LineWidth',0.6,'LineStyle','-'); %#ok<NASGU>

legend_entries = {'Inhomogeneous Def. Shape with NH', ...
    'Posterior Deformed Shape'};

% Plot dummy rectangles for legend
h = zeros(length(legend_entries), 1);
h(1) = plot([0 0], [0 0], 'LineWidth', 0.3, 'LineStyle', '-', 'Color', [0 0.5 0]);
h(2) = plot([0 0], [0 0], 'LineWidth', 0.6, 'LineStyle', '-', 'Color', 'red');
% Show legend
legend(h, legend_entries, 'location', 'northwest', 'Interpreter', 'latex');
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
% Generate PDF and PNG for the paper
axis equal;
box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10)
plot_mean_dispX_Y_contour = ...
    patch(BVP.preProc.msh.xSurf + BVP.statfem.mu_ux_y_pc_Surf, BVP.preProc.msh.ySurf + BVP.statfem.mu_uy_y_pc_Surf, ...
    BVP.statfem.mu_ux_y_pc_Surf, 'FaceColor', 'interp'); %#ok<NASGU>
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
% % axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11)
plot_std_dispX_Y_contour = ...
    patch(BVP.preProc.msh.xSurf + BVP.statfem.mu_ux_y_pc_Surf,BVP.preProc.msh.ySurf+ BVP.statfem.mu_uy_y_pc_Surf, ...
    BVP.statfem.std_ux_y_pc_Surf,'FaceColor','interp'); %#ok<NASGU>
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
% % axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12)
dispError_x = abs(BVP.proc.NH.ux - BVP.statfem.mu_u_y_pc{1});
dispError_x_surf = makeSurf(BVP.preProc.msh.elementNodes,dispError_x);
patch(BVP.preProc.msh.xSurf, BVP.preProc.msh.ySurf, dispError_x_surf, 'FaceColor', 'interp');
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
% % axis off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(13)
dispError_y = abs(BVP.proc.NH.uy - BVP.statfem.mu_u_y_pc{2});
dispError_y_surf = makeSurf(BVP.preProc.msh.elementNodes,dispError_y);
patch(BVP.preProc.msh.xSurf, BVP.preProc.msh.ySurf, dispError_y_surf, 'FaceColor', 'interp');
% Show axis
xlabel('$X$', 'interpreter', 'latex');
ylabel('$Y$', 'interpreter', 'latex');
axis equal; view(0,90); colormap(jet); colorbar('vert');
shading interp;
box on;
% % axis off;
end
