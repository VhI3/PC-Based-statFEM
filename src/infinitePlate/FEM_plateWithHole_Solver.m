function [u, force] = FEM_plateWithHole_Solver(preProc_Variables, E_vector)
% Project: PC-Based-statFEM 
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)

%% Assign from BVP
force                  = preProc_Variables.BC.force;
activeDOFs             = preProc_Variables.BC.activeDOFs;
prescribedDOFs         = preProc_Variables.BC.prescribedDOFs;
u                      = preProc_Variables.proc.LE.u;
GDOFs                  = preProc_Variables.msh.GDOFs;
nElm                   = preProc_Variables.msh.nElm;
nodeCoordinates        = preProc_Variables.msh.nodeCoordinates;
et                     = preProc_Variables.msh.elementNodes;
ed                     = preProc_Variables.msh.elementDOFs;
tickness               = preProc_Variables.geometry.T;
order2D                = preProc_Variables.proc.order2D;
C_Constitutive              = preProc_Variables.material.C_Constitutive;
% CC              = preProc_Variables.material.C_pstress;
gaussLocations         = preProc_Variables.proc.gaussLocations;
gaussWeights           = preProc_Variables.proc.gaussWeights;




%% Calcualtion
% calculation of the system stiffness matrix
stiffness = zeros(GDOFs,GDOFs);% reserve stiffness matrix
for e = 1:nElm   % loop over elements
    indice = et(e,:);
    elementDof = ed(e,:);
    E_e = E_vector(indice);
    ke = zeros(size(elementDof,2));
    for i = 1:size(gaussWeights,1)
        GaussPoint = gaussLocations(i,:);
        xi = GaussPoint(1);
        eta = GaussPoint(2);
        N   = Lagrange2D(xi,eta,1,1,[-1 1],[-1 1]);
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
        E = N*E_e;
        %
        ke = ke + gaussWeights(i)*(B'*E*C_Constitutive*B)*det(JacobianMatrix)*tickness;
    end
    stiffness(elementDof,elementDof)= stiffness(elementDof,elementDof) + ke;
end
% stiffness = globalstiffness_linearElastic2D(preProc_Variables ,E_vector);
%
Krr = sparse(stiffness(activeDOFs,activeDOFs));
Kru = sparse(stiffness(activeDOFs,prescribedDOFs));
Kur = sparse(stiffness(prescribedDOFs,activeDOFs));
Kuu = sparse(stiffness(prescribedDOFs,prescribedDOFs));
%
% Solution
% static analysis
Rr = sparse(force(activeDOFs));
Uu = sparse(u(prescribedDOFs));
Ur = Krr\(Rr-Kru*Uu);
Ru = Kuu*Uu+Kur*Ur;
%
u(activeDOFs) = Ur;
force(prescribedDOFs)= Ru;
%% assign Back to BVP

end
