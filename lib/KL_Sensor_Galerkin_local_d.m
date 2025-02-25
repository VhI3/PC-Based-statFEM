function [eigval,eigvec_f,eigvec,KLterms] = KL_Sensor_Galerkin_local_d(xn_y, partition,Mterms,dom_bound,cov_func,option)
% Project: PC-Based-statFEM 
% Forker: Vahab Narouie
% License: GNU GPL v3.0 (see LICENSE file for details)
% Originally created by Felipe Uribe (Oct/2014)
% Changes and enhancements: Added comments, cleaned code structure,
%                           updated formatting, and improved readability.


%% domain
a = dom_bound{1}(1);
b = dom_bound{1}(2);
%% FEM parameters
ned  = 1;                         % number of dof per node
nen  = 2;                         % number of nodes per element
nnp  = partition;                 % number of nodal points
nfe  = nnp - 1;                   % number of finite elements
neq  = ned*nen;                   % number of element equations
xnod = xn_y';
ndof = ned*nnp;                   % number of degrees of freedom
IEN  = [(1:(nnp-1))' (2:nnp)']';  % loc-glo
ID   = 1:ndof;

% localization matrix
LM = cell(nfe,1);
for e = 1:nfe
    LM{e} = [ID(IEN(1,e)); ID(IEN(2,e))];
end

% two-node 1D elements: lagrangian shape functions
Nshape = @(xi) [ -(xi-1)/2        % N1
    (xi+1)/2 ]';    % N2
EnergyRatio  = 0.99;

%% Gauss-Legendre parameters
nGp         = 2;              % Gauss-Legendre quadrature order
[x_gl,w_gl] = Gauss_int(nGp);
%% Computing matrix B
B = sparse(ndof,ndof);

for e = 1:nfe
    Be     = zeros(neq);
    det_Je = (xnod(IEN(2,e))-xnod(IEN(1,e)))/2;

    % gauss-legendre quadrature integration
    for p = 1:nGp
        xi_gl = x_gl(p);
        NN    = Nshape(xi_gl);    % shape functions on GL points
        % B matrix
        Be = Be + NN'*NN*det_Je*w_gl(p);
    end
    B(LM{e},LM{e}) = B(LM{e},LM{e}) + Be;
end


%% computing matrix C
C = sparse(ndof,ndof);

for e = 1:nfe
    xe     = xnod(IEN(:,e));
    det_Je = (xnod(IEN(2,e))-xnod(IEN(1,e)))/2;
    for f = 1:nfe
        Cef    = zeros(nen);
        xf     = xnod(IEN(:,f));
        det_Jf = (xnod(IEN(2,f))-xnod(IEN(1,f)))/2;
        for p1 = 1:nGp
            xi_gl_e = x_gl(p1);
            NNe     = Nshape(xi_gl_e);   % shape function of element e
            xp1     = sum(NNe.*xe);      % in global coordinates
            for p2 = 1:nGp
                xi_gl_f = x_gl(p2);
                NNf     = Nshape(xi_gl_f);   % shape function of element f
                xp2     = sum(NNf.*xf);       % in global coordinates
                % element C matrix
                %             cov_func(xp1,xp2)
                Cef = Cef + cov_func(xp1,xp2)*NNe'*NNf*det_Je*det_Jf*w_gl(p1)*w_gl(p2);
            end
        end
        C(LM{e},LM{f}) = C(LM{e},LM{f}) + Cef;
    end
end


%% solve generalized eigenvalue problem: CA = BAD
if option == 0
    [D,A]        = eigs(C,B,Mterms);  % or: [D,L] = eig(C,B);
    [eigval,idx] = sort(diag(A),'descend');
    eigvec       = D(:,idx);
elseif option == 1
    %%
    [~,Afull] = eig(full(C),full(B));
    % % % Calculate the cumulated energy ratio
    cumulated_eigs = cumsum(abs(diag(Afull)))/sum(abs(diag(Afull))) ;
    tmp = find(cumulated_eigs >= EnergyRatio) ;
    Mterms = tmp(1) ;
    [D,A]        = eigs(C,B,Mterms);
    [eigval,idx] = sort(diag(A),'descend');
    eigvec       = D(:,idx);
elseif option == 2
    [D,A]        = eigs(C,B,Mterms);
    [eigval,idx] = sort(diag(A),'descend');
    eigvec       = D(:,idx);
end
%%

KLterms = Mterms;
% function form
eigvec_f = cell(KLterms,1);
for k = 1:KLterms
    eigvec_f{k} = @(xx) interp1(xnod,eigvec(:,k),xx);
end

end
