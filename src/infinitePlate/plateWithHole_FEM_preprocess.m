%% plateWithHole_FEM_preprocess.m
% --------------------------------------------------------------
% Preprocess the finite element model for a 2D plate with a hole.
%
% Description:
% This function sets up the geometry, mesh properties, material properties,
% boundary conditions, and numerical integration parameters for a plate with
% a hole under traction loading. It is used for both deterministic and
% stochastic finite element analyses, including uncertainty quantification.
%
% Inputs:
%   BVP : Structure containing mesh information and model parameters.
%
% Outputs:
%   BVP : Updated structure with preprocessed data ready for FEM analysis.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

function BVP = plateWithHole_FEM_preprocess(BVP)
%% Assign from BVP
BVP                         = quadsCheck(BVP);
nodeCoordinates             = BVP.preProc.msh.nodeCoordinates;
elementNodes                = BVP.preProc.msh.elementNodes;
%% Calculation
% -------------------- 1. Geometry properties -----------------------------
% -------------------- 1. Geometry properties -----------------------------
nuElm                       = size(elementNodes,1);
L                           = max(nodeCoordinates(:,1)); % Length of the plate [m]
H                           = max(nodeCoordinates(:,2)); % Height of the plate [m]
R                           = 0.02;      % Radius of the plate
T                           = 0.01;      % Thickness of the plate [m]
DIM                         = 2;         % Dimension of the problem
% -------------------- 2. Mesh properties ---------------------------------
DOFs                       = 2;           % Degree of Freedom per Nodes
nElNode                    = 4;
eDOFs                      = DOFs*nElNode;  % Degree of Freedom per element
nuNodes                    = size(nodeCoordinates,1);
GDOFs                      = nuNodes*DOFs;  % Global Degree of Freedom per Nodes
% Generate every element DOFs
elementDOFs                = zeros(nuElm,eDOFs);
elementDOFs(:,1:2:eDOFs)   = 2*elementNodes-1;
elementDOFs(:,2:2:eDOFs)   = 2*elementNodes;
%
leftNodes                  = find(abs(nodeCoordinates(:,1))< eps);
rightNodes                 = find(abs(nodeCoordinates(:,1)- L)< eps);
topNodes                   = find(abs(nodeCoordinates(:,2)- H)< eps);
bottomNodes                = find(abs(nodeCoordinates(:,2))< eps);
cornerNode                 = intersect(rightNodes,bottomNodes);
cornerNodeDofX             = DOFs*cornerNode-1;
%
bottomNodesDofx1           = DOFs*bottomNodes-1;
bottomNodesDofx            = [bottomNodesDofx1(1,1) ; bottomNodesDofx1(3:end,1); bottomNodesDofx1(2,1)];
bottomNodesCoor1           = nodeCoordinates(bottomNodes,:);
bottomNodesCoor            = [bottomNodesCoor1(1,1) ; bottomNodesCoor1(3:end,1); bottomNodesCoor1(2,1)];
% -------------------- 3. Material properties -----------------------------
% deterministic Young's modulus
E                         = 200e3;
% We added five layer aotund the hole and every layer has its own Young's modulus.
% Young's modulus of the first layer
E_punch1                  = 10e3;
% Young's modulus of the second layer
E_punch2                  = 20e3;
% Young's modulus of the third layer
E_punch3                  = 50e3;
% Young's modulus of the forth layer
E_punch4                  = 80e3;
% Young's modulus of the fifth layer
E_punch5                  = 100e3;
% deterministic Possion's ratio
nu                        = 0.5;
% deterministic Lame's constant 1
lambda                    = @(E) nu*E/(1+nu)/(1-2*nu);
% deterministic Lame's constant 2
mu                        = @(E) E/2/(1+nu);
Dmatrix                   = @(E) [lambda(E)+2*mu(E) lambda(E) 0; lambda(E) lambda(E)+2*mu(E) 0; 0 0 mu(E)];
C_Constitutive =...     % Constitutive matrix for plain stress
  1/(1-nu^2)*[  1      nu          0;
  nu        1          0;
  0        0  (1-nu)/2  ];
% The Young's modulus vector on every node
E_vector                  = E*ones(nuNodes,1);
A10_vector                = E_vector/6;
% The Young's modulus vector on every node included the Young's modulus in the vicinity of punch
E_punch_vector            = E_vector;
% Assign the Young's modulus of every layer.
E_punch_vector(BVP.preProc.msh.nodeFirstLayerElement,1) = E_punch1*ones(size(BVP.preProc.msh.nodeFirstLayerElement,2),1);
E_punch_vector(BVP.preProc.msh.nodeSecondLayerElement,1) = E_punch2*ones(size(BVP.preProc.msh.nodeSecondLayerElement,2),1);
E_punch_vector(BVP.preProc.msh.nodeThirdLayerElement,1) = E_punch3*ones(size(BVP.preProc.msh.nodeThirdLayerElement,2),1);
E_punch_vector(BVP.preProc.msh.nodeFourthLayerElement,1) = E_punch4*ones(size(BVP.preProc.msh.nodeFourthLayerElement,2),1);
E_punch_vector(BVP.preProc.msh.nodeFifthLayerElement,1) = E_punch5*ones(size(BVP.preProc.msh.nodeFifthLayerElement,2),1);

A10_punch_vector = E_punch_vector/6;

% -------------------- 4. Uncertain material properties -------------------
mu_E         = 2e5;                    % mean value of Young's modulus
sig_E        = 3e4;            % standard deviation of Young's modulus
nMCS         = 1e1;         % number of monte carlo simulations once
mu_kappa = log(mu_E ^ 2 / sqrt(mu_E ^ 2 + sig_E ^ 2));
sig_kappa = sqrt(log(1 + sig_E ^ 2 / mu_E ^ 2));
% -------------------- 5. Exact Quantities---------------------------------
%
% -------------------- 6. processing properties ---------------------------
u          = zeros(GDOFs,1);
u_nonHom   = zeros(GDOFs,1);

xiRange1   = [-1 1];    % Can be changed. for example : [0 1]
xiRange2   = [-1 1];    % Can be changed. for example : [0 1]
NGT        = 2;
[xi1, w1]  = Gauss_int(NGT);
[xi2, w2]  = Gauss_int(NGT);
dispSF     = 1;
order2D    = [1 2 4 3];
counter = 1;
for j = 1:NGT
  for i = 1:NGT
    gaussLocations(counter,:)= [xi1(i) xi2(j)];
    gaussWeights(counter,:) = w1(i)*w2(j);
    counter = counter + 1;
  end
end
gaussLocations = gaussLocations(order2D,:);
gaussWeights = gaussWeights(order2D,:);
% -------------------- 7. Boundary conditions -----------------------------
% Dirichlet
prescribedDOFsX = DOFs*leftNodes-1;
prescribedDOFsY = DOFs*leftNodes;

prescribedDOFs = sort([prescribedDOFsX ; prescribedDOFsY]);
% prescriebed displacements
dbc = u(prescribedDOFs);
% number of displacement constrained nodes
ndbc = size(prescribedDOFs,1);
% prescriebed force
traction = 300;
if size(rightNodes,1)>2
  rightNodesCorr = [rightNodes(1); rightNodes(3:end); rightNodes(2)];
else
  rightNodesCorr = rightNodes;
end
rightEdge_y1 = nodeCoordinates(rightNodesCorr,2);
rightEdge_y = [rightEdge_y1(1:end-1) rightEdge_y1(2:end)];
rightEdge_y_dofX1 = DOFs*rightNodesCorr-1;
rightEdge_y_dofX = [rightEdge_y_dofX1(1:end-1) rightEdge_y_dofX1(2:end)];
force = zeros(GDOFs,1);
for ey = 1:size(rightEdge_y,1)
  ed = rightEdge_y_dofX(ey,:); ec = rightEdge_y(ey,:);
  force(ed,1) = force(ed,1) + abs(ec(2)-ec(1))/2*traction*ones(2,1);
end

tbc = [];
tbc = [tbc rightNodes];
ntbc = size(tbc,1);
%
activeDOFsX = setdiff(1:2:GDOFs,prescribedDOFsX);
activeDOFsY = setdiff(2:2:GDOFs,prescribedDOFsY);
activeDOFs  = setdiff(1:GDOFs,prescribedDOFs);
size_activeDOFs = size(activeDOFs,2);
%
P_e = 1;
sig_e = 0.001;
nMC = 1e3;

%% Assign back to BVP
% -------------------- 1. Geometry properties -----------------------------
BVP.preProc.geometry.L               = L;
BVP.preProc.geometry.H               = H;
BVP.preProc.geometry.T               = T;
BVP.preProc.geometry.R               = R;
BVP.preProc.geometry.DIM             = DIM;
% -------------------- 2. Mesh properties ---------------------------------
BVP.preProc.msh.nElm                 = nuElm;
BVP.preProc.msh.nodeCoordinates      = nodeCoordinates;
BVP.preProc.msh.DOFs                 = DOFs;
BVP.preProc.msh.elementNodes         = elementNodes;
BVP.preProc.msh.elementDOFs          = elementDOFs;
BVP.preProc.msh.nuNodes              = nuNodes;
BVP.preProc.msh.eDOFs                = eDOFs;
BVP.preProc.msh.GDOFs                = GDOFs;
BVP.preProc.msh.leftNodes            = leftNodes;
BVP.preProc.msh.rightNodes           = rightNodes;
BVP.preProc.msh.topNodes             = topNodes;
BVP.preProc.msh.bottomNodes          = bottomNodes;
BVP.preProc.msh.bottomNodesDofx      = bottomNodesDofx;
BVP.preProc.msh.bottomNodesCoor      = bottomNodesCoor;
BVP.preProc.msh.cornerNode           = cornerNode;
BVP.preProc.msh.cornerNodeDofX       = cornerNodeDofX;
% -------------------- 3. Material properties -----------------------------
BVP.preProc.material.E               = E;
BVP.preProc.material.E_vector        = E_vector;
BVP.preProc.material.A10_vector      = A10_vector;
BVP.preProc.material.A10_punch_vector= A10_punch_vector;
BVP.preProc.material.nu              = nu;
BVP.preProc.material.lambda           = lambda;
BVP.preProc.material.mu              = mu;
BVP.preProc.material.Dmatrix         = Dmatrix;
BVP.preProc.material.C_Constitutive  = C_Constitutive;
BVP.preProc.material.E_punch1        = E_punch1;
BVP.preProc.material.E_punch2        = E_punch2;
BVP.preProc.material.E_punch3        = E_punch3;
BVP.preProc.material.E_punch4        = E_punch4;
BVP.preProc.material.E_punch5        = E_punch5;
BVP.preProc.material.E_punch_vector  = E_punch_vector;
% -------------------- 4. Uncertain material properties -------------------
BVP.preProc.UQ.mu_E                  = mu_E;
BVP.preProc.UQ.sig_E                 = sig_E;
BVP.preProc.UQ.nMCS                  = nMCS;
BVP.preProc.UQ.mu_kappa              = mu_kappa;
BVP.preProc.UQ.sig_kappa             = sig_kappa;
BVP.preProc.nMC                      = nMC;
% -------------------- 6. processing properties ---------------------------
BVP.preProc.proc.xiRange1            = xiRange1;
BVP.preProc.proc.xiRange2            = xiRange2;
BVP.preProc.proc.order2D             = order2D;
BVP.preProc.proc.dispSF              = dispSF;
BVP.preProc.proc.GP.NGT              = NGT;
BVP.preProc.proc.GP.xi1              = xi1;
BVP.preProc.proc.GP.xi2              = xi2;
BVP.preProc.proc.GP.w11              = w1;
BVP.preProc.proc.GP.w12              = w2;
BVP.preProc.proc.LE.u                = u;
BVP.preProc.proc.LE.u_nonHom         = u_nonHom;
BVP.preProc.proc.gaussWeights        = gaussWeights;
BVP.preProc.proc.gaussLocations      = gaussLocations;
% -------------------- 7. Boundary conditions -----------------------------
BVP.preProc.BC.prescribedDOFs        = prescribedDOFs;
BVP.preProc.BC.prescribedDOFsX       = prescribedDOFsX;
BVP.preProc.BC.prescribedDOFsY       = prescribedDOFsY;
BVP.preProc.BC.activeDOFs            = activeDOFs;
BVP.preProc.BC.activeDOFsX           = activeDOFsX;
BVP.preProc.BC.activeDOFsY           = activeDOFsY;
BVP.preProc.BC.ndbc                  = ndbc;
BVP.preProc.BC.ntbc                  = ntbc;
BVP.preProc.BC.dbc                   = dbc;
BVP.preProc.BC.force                 = force;
BVP.preProc.BC.force_nonHom          = force;
BVP.preProc.BC.size_activeDOFs       = size_activeDOFs;

BVP.preProc.sig_e                    = sig_e;
BVP.preProc.P_e                      = P_e;

disp('1. Pre-processing for the Plate with hole problem is completed.');
end
