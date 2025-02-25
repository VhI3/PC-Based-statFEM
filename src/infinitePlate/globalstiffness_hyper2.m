function [ stiffness, R ] = globalstiffness_hyper2( BVP_tmp, u )
% Project: PC-Based-statFEM 
% Author: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)


%% Assign from BVP
GDOFs = BVP_tmp.GDOFs;
dof   = 2;                 % degree of freedom
nElm  = BVP_tmp.nElm;
nn = 4;  % number of nodes for each element
nodeCoordinates = BVP_tmp.nodeCoordinates;
et    = BVP_tmp.elementNodes;
ed    = BVP_tmp.elementDOFs;
xi1 = BVP_tmp.GP.xi1;
xi2 = BVP_tmp.GP.xi2;
w1  = BVP_tmp.GP.w11;
w2  = BVP_tmp.GP.w12;
nnElem = 8;
A10_vector = BVP_tmp.mat;
thickness = BVP_tmp.thickness;
%% Calculation
% initialazion
stiffness = zeros(GDOFs,GDOFs);             % reserve stiffness matrix
R = zeros(GDOFs,1);               % reserve residual matrix
%
for e = 1:nElm                % loop over elements
    indice = et(e,:);     elementDof = ed(e,:);
    elDisp = u(elementDof);
    elDisp = reshape(elDisp, dof, nn);
    A10_e = A10_vector(indice);
    ke = zeros(size(elementDof,2));
    re = zeros(size(elementDof,2),1);
    for j = 1:size(xi2,2)
        for i = 1:size(xi1,2)
            N    = Lagrange2D(xi1(i),xi2(j),1,1,[-1 1],[-1 1]);
            N1_d = Lagrange2D_d(xi1(i),xi2(j),1,1,[-1 1],[-1 1],1);
            N2_d = Lagrange2D_d(xi1(i),xi2(j),1,1,[-1 1],[-1 1],2);
            N1_d = N1_d([1 2 4 3]);
            N2_d = N2_d([1 2 4 3]);
            naturalDerivatives = [N1_d;N2_d]';
            %
            [JacobianMatrix,~,XYDerivatives] = Jacobian2D(nodeCoordinates(indice,:),naturalDerivatives);
            F = elDisp * XYDerivatives + eye(dof);
            
            C_tensor = F'*F;      % right Cauchy-Green deformation tensor
            C11 = C_tensor(1,1); C12 = C_tensor(1,2); C21 = C_tensor(2,1); C22 = C_tensor(2,2);
            C33 = 1 / ( C11*C22 - C12*C21 );   % satisfy det(C_tensor) = 1 for incompressible material with Poisson ratio = 0.5
            I1C = [1-C33^2*C22; 1-C33^2*C11; C33^2*C12];
            I1CC = C33^2*[     2*C33*C22^2,  2*C33*C11*C22-1,    -2*C33*C22*C12;
                2*C33*C11*C22-1,      2*C33*C11^2,    -2*C33*C11*C12;
                -2*C33*C22*C12,   -2*C33*C11*C12,   2*C33*C12^2+0.5;];
            
            A10 = N*A10_e;
            
            stress = 2*A10*I1C;
            dtan = 4*A10*I1CC;
            
            %
            BN = zeros(3,nnElem);
            BG = zeros(4,nnElem);
            for k = 1:nn
                BN(:,k*2-1:k*2) = [ F(1,1)*XYDerivatives(k,1)     F(2,1)*XYDerivatives(k,1);
                    F(1,2)*XYDerivatives(k,2)     F(2,2)*XYDerivatives(k,2);
                    F(1,1)*XYDerivatives(k,2)+ F(1,2)*XYDerivatives(k,2)  F(2,1)*XYDerivatives(k,2)+F(2,2)*XYDerivatives(k,1)];
                
                BG(:,k*2-1:k*2) = [XYDerivatives(k,1)          0;
                    XYDerivatives(k,2)          0;
                    0   XYDerivatives(k,1);
                    0   XYDerivatives(k,2);];
                
            end
            sigma = [stress(1) stress(3);
                stress(3) stress(2)];
            stan = zeros(4);
            stan(1:2,1:2) = sigma;
            stan(3:4,3:4) = sigma;
            %
            ke = ke + w1(i)*w2(j)*(BN'*dtan*BN + BG'*stan*BG)*det(JacobianMatrix)*thickness;
            re = re + w1(i)*w2(j)*BN'*stress*det(JacobianMatrix)*thickness;
        end
    end
    stiffness(elementDof,elementDof)= stiffness(elementDof,elementDof) + ke;
    R(elementDof) = R(elementDof) + re;
end
end
