function BVP = quadsCheck(BVP)
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)


%% Assig from BVP
% Add the generated mesh from gmsh to Matlab
plate_with_hole_Mesh % 2D algorithm = Delaunay, 3D algorithm = Frontal, 2D recombination algorithm = Blossom, Recombine all trangular meshes, Subdivision algorithm = All Quads, Element size factor = 10

%% Calculation
% Extract node coordinates excluding the last column which might be z-coordinates or other data
nodeCoordinates = msh.POS(:, 1:end-1);

% Extract element node indices excluding the last column which might be additional data
elementNodes = msh.QUADS(:, 1:end-1);

% Define orders for Jacobian calculation and for 2D shape functions
% This because of the lagrange shape functions I programmed in my Master.
order_jacobi = [1 4 3 2];
order2D = [1 2 4 3];

% Define natural derivatives of shape functions in reference space
N1_d = [-0.25, 0.25, -0.25, 0.25];
N2_d = [-0.25, -0.25, 0.25, 0.25];

% Reorder the natural derivatives according to the specified order for 2D elements
N1_d = N1_d(order2D);
N2_d = N2_d(order2D);

% Combine the natural derivatives into a single matrix
naturalDerivatives = [N1_d; N2_d]';

% Initialize a counter for elements with negative Jacobian determinants
count = 0;

% -------------------------------------------------------------------------
% Element Layer Classification Around the Hole (Young's Modulus Distribution)
% -------------------------------------------------------------------------
% The plate contains **five concentric layers** surrounding the hole, each
% assigned a distinct Young's modulus. Elements are classified into layers
% based on their area (totalArea), representing varying material stiffness.
%
% Young's Modulus Assignments:
%   Layer 1 (Innermost):   E =  20 GPa → totalArea < 0.25e-4
%   Layer 2:              E =  40 GPa → 0.25e-4 ≤ totalArea < 0.5e-4
%   Layer 3:              E =  80 GPa → 0.5e-4  ≤ totalArea < 0.6e-4
%   Layer 4:              E = 120 GPa → 0.6e-4  ≤ totalArea < 0.8e-4
%   Layer 5:              E = 160 GPa → 0.8e-4  ≤ totalArea < 1.5e-4
%   Layer 6 (Outer):      E = 200 GPa → totalArea ≥ 1.5e-4
%
% --------------------------
% Purpose of the Layering:
% - Simulates material heterogeneity around the hole.
% - Captures the effect of varying stiffness on stress and displacement fields.
% - Enhances model fidelity for statFEM analysis.
%
% --------------------------
% How Classification Works:
% - Each quadrilateral element’s area is calculated by dividing it into two triangles.
% - The computed `totalArea` determines which layer the element belongs to.
% - Elements are stored in separate arrays corresponding to their layers.
% -------------------------------------------------------------------------

%
% This variable is used to store the indices of small elements
firstLayerElementIndices = [];
firstLayerElementNodes = [];
secondLayerElementIndices = [];
secondLayerElementNodes = [];
thirdLayerElementIndices = [];
thirdLayerElementNodes = [];
fourthLayerElementIndices = [];
fourthLayerElementNodes = [];
fifthLayerElementIndices = [];
fifthLayerElementNodes = [];
sixthLayerElementIndices = [];
sixthLayerElementNodes = [];

smallElementIndices = []; %#ok<NASGU> % Initialize an empty array to store indices of small elements
smallElementNodes = []; %#ok<NASGU> % Initialize an empty array to store nodes of small elements

% Loop through each element in the mesh
for e = 1:size(elementNodes, 1)
  % Extract node indices for the current element
  indice = elementNodes(e, :);
  
  % Get the coordinates of the nodes for the current element
  coords = nodeCoordinates(elementNodes(e, :), :);
  
  % Calculate the Jacobian matrix using the node coordinates and natural derivatives
  [JacobianMatrix, ~, ~] = Jacobian2D(coords, naturalDerivatives);
  
  % Calculate the determinant of the Jacobian matrix
  detJ = det(JacobianMatrix);
  
  % If the determinant is negative, adjust the node ordering
  if detJ < 0
    count = count + 1;  % Increment the counter
    elementNodes(e, :) = indice(order_jacobi);  % Reorder the nodes
  end
  
  % For a quadrilateral, calculate the area by dividing it into two triangles
  % Triangle 1: nodes 1, 2, 3
  % Triangle 2: nodes 1, 3, 4
  % You can use the Shoelace formula or any other method you prefer
  area1 = polyarea(coords([1, 2, 3], 1), coords([1, 2, 3], 2));
  area2 = polyarea(coords([1, 3, 4], 1), coords([1, 3, 4], 2));
  
  % Sum the areas of the two triangles to get the area of the quadrilateral
  totalArea = area1 + area2;
  
  % Assign this element to the appropriate layer
  if totalArea < 0.25e-4
    firstLayerElementIndices = [firstLayerElementIndices, e]; %#ok<AGROW> % Add element index to the list
    firstLayerElementNodes = [firstLayerElementNodes; indice']; %#ok<AGROW>
  elseif (totalArea < 0.5e-4 && totalArea > 0.25e-4)
    secondLayerElementIndices = [secondLayerElementIndices, e]; % Add element index to the list
    secondLayerElementNodes = [secondLayerElementNodes; indice'];
  elseif (totalArea < 0.6e-4 && totalArea > 0.5e-4)
    thirdLayerElementIndices = [thirdLayerElementIndices, e]; % Add element index to the list
    thirdLayerElementNodes = [thirdLayerElementNodes; indice'];
  elseif (totalArea < 0.8e-4 && totalArea > 0.6e-4)
    fourthLayerElementIndices = [fourthLayerElementIndices, e]; % Add element index to the list
    fourthLayerElementNodes = [fourthLayerElementNodes; indice'];
  elseif (totalArea < 1.5e-4 && totalArea > 0.8e-4)
    fifthLayerElementIndices = [fifthLayerElementIndices, e]; % Add element index to the list
    fifthLayerElementNodes = [fifthLayerElementNodes; indice'];
  else
    sixthLayerElementIndices = [sixthLayerElementIndices, e]; % Add element index to the list
    sixthLayerElementNodes = [sixthLayerElementNodes; indice'];
  end
end

% nodeSmallElement = unique(smallElementNodes);
nodeFirstLayerElement = unique(firstLayerElementNodes);
nodeSecondLayerElement = unique(secondLayerElementNodes);
nodeThirdLayerElement = unique(thirdLayerElementNodes);
nodeFourthLayerElement = unique(fourthLayerElementNodes);
nodeFifthLayerElement = unique(fifthLayerElementNodes);

% Create a subplot which plots the mesh, and sensor locations
xSurf = makeSurf(elementNodes,nodeCoordinates(:,1));
ySurf = makeSurf(elementNodes,nodeCoordinates(:,2));

% Importing the sensor data, choose based on the number of sensors
% sensor_nsen112 % nSen = 112 % Paper relevant
sensor_nsen32 % nSen = 32 % Paper relevant
% sensor_nsen11 % nSen = 11 % Paper relevant

% Extract sensor coordinates excluding the last column which might be z-coordinates or other data
%
sensorCoordinates        = [sensor.POS(:, 1) sensor.POS(:, 2)];
[sensorIndex,~]          = dsearchn(nodeCoordinates,sensorCoordinates);
sensorCoordinates        = [nodeCoordinates(sensorIndex,1) nodeCoordinates(sensorIndex,2)];
% This part is to patch the six layers in the mesh.
xSurf_6 = makeSurf(elementNodes(sixthLayerElementIndices,:),nodeCoordinates(:,1));
ySurf_6 = makeSurf(elementNodes(sixthLayerElementIndices,:),nodeCoordinates(:,2));
%
xSurf_5 = makeSurf(elementNodes(fifthLayerElementIndices,:),nodeCoordinates(:,1));
ySurf_5 = makeSurf(elementNodes(fifthLayerElementIndices,:),nodeCoordinates(:,2));
% hold on
% patch(xSurf_5,ySurf_5,'black','FaceAlpha',.1,'EdgeColor','black','LineWidth',.2,'LineStyle',':');
xSurf_4 = makeSurf(elementNodes(fourthLayerElementIndices,:),nodeCoordinates(:,1));
ySurf_4 = makeSurf(elementNodes(fourthLayerElementIndices,:),nodeCoordinates(:,2));
%
xSurf_3 = makeSurf(elementNodes(thirdLayerElementIndices,:),nodeCoordinates(:,1));
ySurf_3 = makeSurf(elementNodes(thirdLayerElementIndices,:),nodeCoordinates(:,2));
%
xSurf_2 = makeSurf(elementNodes(secondLayerElementIndices,:),nodeCoordinates(:,1));
ySurf_2 = makeSurf(elementNodes(secondLayerElementIndices,:),nodeCoordinates(:,2));
%
xSurf_1 = makeSurf(elementNodes(firstLayerElementIndices,:),nodeCoordinates(:,1));
ySurf_1 = makeSurf(elementNodes(firstLayerElementIndices,:),nodeCoordinates(:,2));

sensorNodes = sensor.QUADS(:, 1:end-1);
nuSenElm                       = size(sensorNodes,1);
nuSensors                    = size(sensorCoordinates,1);
% % Plot the sensors as elements
x_Surf_sensors = makeSurf(sensorNodes,sensorCoordinates(:,1));
y_Surf_sensors = makeSurf(sensorNodes,sensorCoordinates(:,2));

% Find the sensors on the dirchlet boundary
sensorNull = find(abs(sensorCoordinates(:,1))< eps);


%% Assign back to BVP
BVP.preProc.msh.nodeCoordinates           = nodeCoordinates;
BVP.preProc.msh.elementNodes              = elementNodes;
BVP.preProc.msh.xSurf                     = xSurf;
BVP.preProc.msh.ySurf                     = ySurf;
BVP.preProc.msh.firstLayerElementIndices  = firstLayerElementIndices;
BVP.preProc.msh.secondLayerElementIndices = secondLayerElementIndices;
BVP.preProc.msh.thirdLayerElementIndices  = thirdLayerElementIndices;
BVP.preProc.msh.fourthLayerElementIndices = fourthLayerElementIndices;
BVP.preProc.msh.fifthLayerElementIndices  = fifthLayerElementIndices;
BVP.preProc.msh.xSurf_1                   = xSurf_1;
BVP.preProc.msh.ySurf_1                   = ySurf_1;
BVP.preProc.msh.xSurf_2                   = xSurf_2;
BVP.preProc.msh.ySurf_2                   = ySurf_2;
BVP.preProc.msh.xSurf_3                   = xSurf_3;
BVP.preProc.msh.ySurf_3                   = ySurf_3;
BVP.preProc.msh.xSurf_4                   = xSurf_4;
BVP.preProc.msh.ySurf_4                   = ySurf_4;
BVP.preProc.msh.xSurf_5                   = xSurf_5;
BVP.preProc.msh.ySurf_5                   = ySurf_5;
BVP.preProc.msh.xSurf_6                   = xSurf_6;
BVP.preProc.msh.ySurf_6                   = ySurf_6;
BVP.preProc.msh.nodeFirstLayerElement     = nodeFirstLayerElement;
BVP.preProc.msh.nodeSecondLayerElement    = nodeSecondLayerElement;
BVP.preProc.msh.nodeThirdLayerElement     = nodeThirdLayerElement;
BVP.preProc.msh.nodeFourthLayerElement    = nodeFourthLayerElement;
BVP.preProc.msh.nodeFifthLayerElement     = nodeFifthLayerElement;


% BVP.preProc.msh.xSurf_smallElement    = xSurf_smallElement;
% BVP.preProc.msh.ySurf_smallElement    = ySurf_smallElement;
% BVP.preProc.msh.smallElementIndices   = smallElementIndices;
% BVP.preProc.msh.nodeSmallElement      = nodeSmallElement;

BVP.preProc.statfem.sensorCoordinates = sensorCoordinates;
BVP.preProc.statfem.sensorNodes       = sensorNodes;
BVP.preProc.statfem.nuSenElm          = nuSenElm;
BVP.preProc.statfem.nuSensors         = nuSensors;
BVP.preProc.statfem.sensorNull        = sensorNull;
BVP.preProc.statfem.x_Surf_sensors = x_Surf_sensors;
BVP.preProc.statfem.y_Surf_sensors = y_Surf_sensors;

end
