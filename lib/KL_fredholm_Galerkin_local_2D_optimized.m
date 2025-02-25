function [eigval, eigvec_f, eigvec, KLterms] = KL_fredholm_Galerkin_local_2D_optimized(preProc_Variables, cov_func)
% Project: PC-Based-statFEM
% Forker: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% Originally created by Felipe Uribe (Oct/2014)
% Changes and enhancements: Added comments, cleaned code structure,
%                           updated formatting, and improved readability.
% --------------------------------------------------------------

% Unpack variables from preProc_Variables
numberElements  = preProc_Variables.msh.nElm;
nodeCoordinates = preProc_Variables.msh.nodeCoordinates;
et              = preProc_Variables.msh.elementNodes;
% nuNodes         = preProc_Variables.msh.nuNodes;
order2D         = preProc_Variables.proc.order2D;

% Initialization of shape functions and their derivatives
Nshape = @(xi, eta) [(eta - 1)*(xi - 1), -(eta - 1)*(xi + 1), (eta + 1)*(xi + 1), -(eta + 1)*(xi - 1)] / 4;
dN_dxi = @(xi, eta) [(eta - 1), -(eta - 1), (eta + 1), -(eta + 1)] / 4;
dN_deta = @(xi, eta) [(xi - 1), -(xi + 1), (xi + 1), -(xi - 1)] / 4;

% Setting up Gauss-Legendre quadrature
nGp = 2;
[xi1, w1] = Gauss_int(nGp);
[xi2, w2] = Gauss_int(nGp);
gaussLocations = zeros(nGp^2, 2);
gaussWeights = zeros(nGp^2, 1);

% Compute the tensor product of the 1D Gauss points and weights
counter = 1;
for j = 1:nGp
    for i = 1:nGp
        gaussLocations(counter, :) = [xi1(i), xi2(j)];
        gaussWeights(counter) = w1(i) * w2(j);
        counter = counter + 1;
    end
end
gaussLocations = gaussLocations(order2D,:);
gaussWeights = gaussWeights(order2D,:);

% Optimized computation for matrices B and C
[B, C] = computeMatricesBC(numberElements, et, nodeCoordinates, Nshape, dN_dxi, dN_deta, gaussLocations, gaussWeights, cov_func);

% Solve the generalized eigenvalue problem: CA = BAD
[~, Afull] = eig(C, B);
cumulated_eigs = cumsum(abs(diag(Afull))) / sum(abs(diag(Afull)));
KLterms = find(cumulated_eigs >= 0.99, 1);
[D, A] = eigs(sparse(C), sparse(B), KLterms);
[eigval, idx] = sort(diag(A), 'descend');
eigvec = D(:, idx);

% Constructing function handles for eigenfunctions
eigvec_f = arrayfun(@(k) @(xx) interp1(nodeCoordinates, eigvec(:, k), xx), 1:KLterms, 'UniformOutput', false);
end

function [B, C] = computeMatricesBC(numberElements, et, nodeCoordinates, Nshape, dN_dxi, dN_deta, gaussLocations, gaussWeights, cov_func)
nuNodes = size(nodeCoordinates, 1);
gaussWeightsSize = size(gaussWeights, 1);

% Precompute the shape functions and their derivatives for all Gauss points
N_all = arrayfun(@(i) Nshape(gaussLocations(i, 1), gaussLocations(i, 2)), 1:gaussWeightsSize, 'UniformOutput', false);
N_dxi_all = arrayfun(@(i) dN_dxi(gaussLocations(i, 1), gaussLocations(i, 2)), 1:gaussWeightsSize, 'UniformOutput', false);
N_deta_all = arrayfun(@(i) dN_deta(gaussLocations(i, 1), gaussLocations(i, 2)), 1:gaussWeightsSize, 'UniformOutput', false);

B = zeros(nuNodes, nuNodes);
C = zeros(nuNodes, nuNodes);

% Assemble matrix B
for e = 1:numberElements
    indice = et(e, :);
    nodeCoords = nodeCoordinates(indice, :);
    Be = zeros(length(indice));

    for i = 1:gaussWeightsSize
        N = N_all{i};
        N_dxi = N_dxi_all{i};
        N_deta = N_deta_all{i};
        naturalDerivatives = [N_dxi; N_deta]';
        [JacobianMatrix, ~, ~] = Jacobian2D(nodeCoords, naturalDerivatives);

        Be = Be + gaussWeights(i) * (N' * N) * det(JacobianMatrix);
    end
    B(indice, indice) = B(indice, indice) + Be;
end

% Assemble matrix C
for e = 1:numberElements
    indice_e = et(e, :);
    nodeCoords_e = nodeCoordinates(indice_e, :);

    for f = 1:numberElements
        indice_f = et(f, :);
        nodeCoords_f = nodeCoordinates(indice_f, :);
        Cef = zeros(length(indice_f));

        for i = 1:gaussWeightsSize
            Ne = N_all{i};
            N1e_d = N_dxi_all{i};
            N2e_d = N_deta_all{i};
            naturalDerivatives_e = [N1e_d; N2e_d]';
            [JacobianMatrix_e, ~, ~] = Jacobian2D(nodeCoords_e, naturalDerivatives_e);
            xp1 = Ne * nodeCoords_e;

            for j = 1:gaussWeightsSize
                Nf = N_all{j}; % Since Nf depends on j, not precomputed for f
                naturalDerivatives_f = [N_dxi_all{j}; N_deta_all{j}]';
                [JacobianMatrix_f, ~, ~] = Jacobian2D(nodeCoords_f, naturalDerivatives_f);
                xp2 = Nf * nodeCoords_f;

                Cef = Cef + gaussWeights(i) * cov_func(xp1, xp2) * (Ne' * Nf) * det(JacobianMatrix_e) * det(JacobianMatrix_f);
            end
        end
        C(indice_e, indice_f) = C(indice_e, indice_f) + Cef;
    end
end
end
