%% plateWithHole_FEM_processNH.m
% --------------------------------------------------------------
% Nonlinear Hyperelastic (Neo-Hookean) FEM Analysis for a Plate with a Hole.
%
% Description:
% This function performs a nonlinear FEM analysis of a plate with a hole 
% using a Neo-Hookean material model. The solution is computed via an 
% incremental-iterative procedure with convergence checks.
%
% Inputs:
%   BVP : Structure containing mesh, material properties, and boundary conditions.
%
% Outputs:
%   BVP : Updated structure with nodal displacements, iteration history, and surface fields.
%
% Project: PC-Based-statFEM
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0 (see LICENSE file for details)
% --------------------------------------------------------------

function BVP = plateWithHole_FEM_processNH(BVP)
%% Assign from BVP
GDOFs           = BVP.preProc.msh.GDOFs;

tol                          = 1e-10;    % tolerance of convergence
maxit                        = 40;      % maximum iterative steps
reit                         = 0;       % reduction index
maxreit                      = 6;% maximum load/displacement step reduction
% timeInterval                 = 0.5;  % initial time interval
timeInterval                 = 1;  % initial time interval


ndbc            = BVP.preProc.BC.ndbc;
ntbc            = BVP.preProc.BC.ntbc;
dbc             = BVP.preProc.BC.dbc;
prescribedDOFs  = BVP.preProc.BC.prescribedDOFs;
force           = BVP.preProc.BC.force;

%
BVP_tmp.GP              = BVP.preProc.proc.GP;
BVP_tmp.nElm            = BVP.preProc.msh.nElm;
BVP_tmp.nodeCoordinates = BVP.preProc.msh.nodeCoordinates;
BVP_tmp.elementNodes    = BVP.preProc.msh.elementNodes;
BVP_tmp.elementDOFs     = BVP.preProc.msh.elementDOFs;
BVP_tmp.GDOFs           = BVP.preProc.msh.GDOFs;
BVP_tmp.thickness        = BVP.preProc.geometry.T;
BVP_tmp.mat             = BVP.preProc.material.A10_punch_vector;
%% Calculation
u = zeros(GDOFs,1);   % nodal displacements
u_steps = u;
cu = zeros(GDOFs,1);  % converged nodal displacements
step = 0;             % load/displacement step index
curtime = 0;          % current time
%
%
cnit = [];            % record the iterative steps
counter = 0;
while curtime ~= 1    % get to the end
  counter = counter + 1;
  curtime = curtime + timeInterval;
  if curtime > 1
    timeInterval = 1 - curtime + timeInterval;
    curtime = 1;
  end
  err = 1e6;        % record iterative error
  perr = err;
  %     num  = err;
  nit = 0;          % record iterative times
  
  iu = zeros(GDOFs,1);   % record iterative displacements
  while (err > tol) && ( nit <= maxit)
    nit = nit+1;
    [ k, r ] = globalstiffness_hyper2( BVP_tmp, u );
    f = zeros(GDOFs,1);         % define external force
    if ntbc~=0                  % enforce traction conditions
      f = force;
    end
    
    if ndbc~=0                  % enforce displacement conditions
      k(prescribedDOFs,:) = zeros(ndbc, GDOFs);
      k(prescribedDOFs,prescribedDOFs) = eye(ndbc);
      f(prescribedDOFs,:) = 0;
      if nit == 1
        f(prescribedDOFs,:) = dbc;
      end
    end
    b = curtime*f - r;          % define right side of the governing equation
    if ndbc~=0
      b(prescribedDOFs) = curtime*dbc - u(prescribedDOFs);
    end
    du = k\b;                   % solve equation
    %
    alldof = 1:GDOFs;
    freedof = setdiff(alldof, prescribedDOFs);    % nodes without displacement constraint
    u = u + du;                 % update displacement
    iu = iu + du;               % update increment displacement
    if nit > 1                  % compute iterative error
      num = b(freedof)' * b(freedof);
      denom = 1+f(freedof)' * f(freedof);
      err = num/denom;
    end
    % output current time step and iterative error
    if err/perr > 1E6 && nit > 2
      nit = maxit+1;   % lf solution diverge extremely, go to next iteration
    else
      perr = err;
    end
  end
  % if convergence
  if  nit <= maxit               % converged
    reit = 0;                  % reset reduction index
    step = step + 1;           % increase converged steps by 1
    cu = u;                    % update converged displacement
    cnit = [cnit, nit]; %#ok<AGROW>
    if length(cnit) >=2 && all(cnit(end-1:end) <= 5)
      timeInterval = timeInterval*1.5;  % increase the increment by times 1.5
    end
    %% save the converged displacement
    u_steps = [u_steps u]; %#ok<AGROW>
    
  else                           % not converged
    if reit <= maxreit         % refine time interval and continue iterating
      curtime = curtime - timeInterval;   % recover current time step
      timeInterval = timeInterval/4;      % refine time interval
      reit = reit+1;         % increase reduction index
      u = cu;                % recover current displacement from last converged displacement
    else
      return;                % stop analysis
    end
  end
end
%
ux_steps                = u_steps(1:2:end,:);
uy_steps                = u_steps(2:2:end,:);
u                       = u_steps(:,end);
ux                      = u(1:2:end);
uy                      = u(2:2:end);
ux_Surf                 = makeSurf(BVP_tmp.elementNodes,ux);
uy_Surf                 = makeSurf(BVP_tmp.elementNodes,uy);

%% Assign back to BVP
BVP.proc.NH.u_steps                   = u_steps;
BVP.proc.NH.ux_steps                  = ux_steps;
BVP.proc.NH.uy_steps                  = uy_steps;
BVP.proc.NH.u                         = u;
BVP.proc.NH.ux                        = ux;
BVP.proc.NH.uy                        = uy;
BVP.proc.NH.ux_Surf                   = ux_Surf;
BVP.proc.NH.uy_Surf                   = uy_Surf;
BVP.proc.ST.b                         = b;


disp('3. Neo-Hooken FEM analysis for the Plate with a Hole is completed.');
end
