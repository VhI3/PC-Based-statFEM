function [ stiffness ] = globalstiffness_linearElastic2D( preproc_variables, E_vector )
% Project: PC-Based-statFEM 
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)


%% Assign from BVP
GDOFs       = preproc_variables.msh.GDOFs;
nElm        = preproc_variables.msh.nElm;
nodeCoordinates = preproc_variables.msh.nodeCoordinates;
et    = preproc_variables.msh.elementNodes;
ed    = preproc_variables.msh.elementDOFs;
tickness = preproc_variables.geometry.T;
order2D = preproc_variables.proc.order2D;
C_pstrain   = preproc_variables.material.C_pstrain;
gaussLocations = preproc_variables.proc.gaussLocations;
gaussWeights   = preproc_variables.proc.gaussWeights;

%% Calculation
stiffness = zeros(GDOFs,GDOFs);% reserve stiffness matrix
for e = 1:nElm   % loop over elements
    indice = et(e,:);
    elementDof = ed(e,:);
    E = E_vector(e);
    ke = zeros(size(elementDof,2));
    for i = 1:size(gaussWeights,1)
        GaussPoint = gaussLocations(i,:);
        xi = GaussPoint(1);
        eta = GaussPoint(2);
        N1_d = Lagrange2D_d(xi,eta,1,1,[-1 1],[-1 1],1);
        N2_d = Lagrange2D_d(xi,eta,1,1,[-1 1],[-1 1],2);
        N1_d = N1_d(order2D);
        N2_d = N2_d(order2D);
        naturalDerivatives = [N1_d;N2_d]';
        %
        [JacobianMatrix,~,XYDerivatives] = Jacobian2D(nodeCoordinates(indice,:),naturalDerivatives);
        % B matrix
        B = zeros(3,size(elementDof,2));
        B(1,1:2:end)    = XYDerivatives(:,1)';
        B(2,2:2:end)    = XYDerivatives(:,2)';
        B(3,1:2:end)    = XYDerivatives(:,2)';
        B(3,2:2:end)    = XYDerivatives(:,1)';
        %
        ke = ke + gaussWeights(i)*(B'*E*C_pstrain*B)*det(JacobianMatrix)*tickness;
    end
    stiffness(elementDof,elementDof)= stiffness(elementDof,elementDof) + ke;
end
end
