%  VNss: Function to calculate quadratic cost
%
%  V_N(theta) = 1/N\sum_{t=1}^N ||y_t - y^p_t||^2
%
%  associated with the bilinear state-space model
%
%  x_{t+1} = Ax_t + F(u_t \kron x_t) + Bu_t + Ke_t
%      y_t = Cx_t + G(u_t \kron x_t) + Du_t + e_t
%
%  where e_t is white noise.  The derivates and Hessian of this cost V_N
%  with respect to a data driven local co-ordinate parametrization theta
%  are also computed.  This function is only used for the estimation of
%  state space model structures.
%
%  Usage is
%
%  [cost,pe,grad,hess] = VNss(Z,theta,OPT,M,div);
%
%  where
%
%   Z         = Input output data in the form specified in introduction
%               to est.m, Type "help est" for details.
%  theta      = Specification of point with respect to vector theta
%               at thich to calculate cost, gradient and Hessian.
%  M          = Data structure which defines the model structure.  For
%               the purpose of this function, which only deals with state
%               space model structures, the elements M.ss.A,--M.ss.K are
%               the only relevant ones.  The are matrices that define the
%               state space model structure given above.
%  OPT        = Data structure which defines options for the estimation
%               algorithm as follows.  It is not used in this function,
%               but is kept as an argument to make it compatible with the
%               general damped Gauss-Newton line search algorithm argmin.m
%  div        = flag variable, that if = 1 causes gradients and hessians
%               to be calculated, but not otherwise; in this latter case
%               only a cost is computed.
%
%    cost     = value of quadratic cost V_N(M).
%    pe       = prediction error sequence.
%    grad     = gradient of V_N(M) with respect to theta
%    hess     = estimate of Hessian of V_N(M) with respect to theta.
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%        		              Australia.

%   Copyright (C) Brett Ninness, Adrian Wills

function [cost,pe,grad,phi,map] = VNss(Z,theta,OPT,M,div)

% Extract input and output from data matrix
[y,u,ny,nu,N] = Z2data(Z);
m=nu; 

% Include delays specified in model structure on inputs
for r=1:nu,
	u(:,r) = [zeros(M.delay(r),1);u(1:N-M.delay(r),r)];
end

% Put input into the correct form
if isempty(u) || nu==0,
 u = zeros(N,0);
end

% Find out if there is a structure that we must obey
if isfield(M,'theta_struct'),
 theta_struct = M.theta_struct;
else
 theta_struct = find(ones(length(theta),1));
end

% Use theta to figure out state space model
M  = theta2m(theta,M);
Ms = M; %save this M for the sampling case

% If continuous-time model required then compute model suitable
% for simulating that system
if M.op=='s',
 if ~isfield(M,'sample'),
  M.sample = 'samplek';
 end
 % Get the model suitable for simulation 
 M = feval(M.sample,M);
end

% Extract system matrices and sizes
A = M.ss.A;
n = size(A,1);
C = M.ss.C;
p = size(C,1);
if isfield(M.ss,'B'),  B  = M.ss.B;  else B  = []; end
if isfield(M.ss,'D'),  D  = M.ss.D;  else D  = []; end
if isfield(M.ss,'K'),  K  = M.ss.K;  else K  = []; end
if isfield(M.ss,'F'),  F  = M.ss.F;  else F  = []; end
if isfield(M.ss,'G'),  G  = M.ss.G;  else G  = []; end
if isfield(M.ss,'X1'), X1 = M.ss.X1; else X1 = []; end

% Determine if any matrices are empty
if numel(B)==0,  B  = zeros(n,nu);   isB  = 0; else isB  = 1; end
if numel(D)==0,  D  = zeros(p,nu);   isD  = 0; else isD  = 1; end
if numel(K)==0,  K  = zeros(n,p);    isK  = 0; else isK  = 1; end
if numel(F)==0,  F  = zeros(n,n*nu); isF  = 0; else isF  = 1; end
if numel(G)==0,  G  = zeros(p,n*nu); isG  = 0; else isG  = 1; end
if numel(X1)==0, X1 = zeros(n,1);    isX1 = 0; else isX1 = 1; end

% Make a small correction if isF and isG are both zero (i.e. linear)
if ~isF && ~isG,
 M.type='ss';
end

% Simulate the state record
if strcmpi(M.type,'ss'),  % Case of linear state space system
 uy  = [u y];
 Ct  = C.';
 Dt  = D.';
 xh  = ltitr(A-K*C,[B-K*D K],uy,X1);
 yh  = xh*Ct + u*Dt;
 pe  = y-yh;
 if div, uyx = [uy xh]; end
else
 % Case of bilinear system - can't just use ltitr to forward simulate that
 [xh,yh,pe,ukx] = systr(y,u,A,B,C,D,K,F,G,X1,isD,isK,isF,isG,isX1,m);
 if div, bigu = [u xh ukx pe]; end
end

% Compute cost
Lam  = (pe'*pe)/N;
switch OPT.cost,
 case 'det'
  [sqrtLam,ispos] = chol(Lam);
  if ~ispos,
   cost=real(prod(diag(sqrtLam))^2);
   if abs(cost)<eps, cost=inf; end
   if isnan(cost), cost=inf; end
   if cost<0, cost=inf; end
  else
   cost = inf;
  end
  if abs(cost)~=inf,
   sqrtLam = sqrtLam\eye(size(Lam));
   scal    = sqrt(cost);
   pe      = pe*sqrtLam*scal;
  end
 otherwise
  cost=trace(Lam);
end

% Compute derivatives by calling VNss_sub() - computes Jacobian matrix of pe
if div,
 deps = eps^(2/3);
 dy   = zeros(1,p);
 dx   = zeros(1,n);
 dxo  = zeros(1,n);
 
 if strcmpi(M.par,'ddlc'),
  Md  = ddlc(theta,Ms);
  nQp = size(Md.Qp,2);
  %map = eye(length(theta_struct));
  map = Md.Qp;
  if isX1, 
   nth = nQp + n; 
   map = blkdiag(map,eye(n));
  else
   nth = nQp; 
  end
 else
  nth = length(theta_struct);
 end
 phi  = zeros(p*N,nth);
 
 % Reserve RAM to store some matrices
 dA = zeros(n);
 dB = zeros(n,nu);
 dC = zeros(p,n);
 dD = zeros(p,nu);
 dK = zeros(n,p);
 dF = zeros(n,n*nu);
 dG = zeros(p,n*nu);
 dx = zeros(n,1);
 
 % Loop over the parameters we are allowed to vary
 for i=1:nth,
  % Compute derivative of system matrices wrt theta based on
  % parametrization (note that op='s' or M.par='grey' will result in numerical
  % derivatives being used).
  if strcmpi(M.op,'s') || strcmpi(M.par,'grey'),
   % Use finite difference to get derivative estimate of state-space
   % model w.r.t. theta (this will often be exact to machine eps)
   thn    = zeros(size(theta));
   thn(theta_struct(i)) = deps;
   switch lower(OPT.sysnd),
    case 'forward' %Use forward numerical differentiation
     pM = theta2m(theta+thn,Ms);
     % If continuous model required then sample these systems
     if M.op=='s', pM  = samplek(pM); end
     nM = M;
     sc = 1;
    case 'backward' % Use backward numerical differentiation
     nM = theta2m(theta-thn,Ms);
     % If continuous model required then sample these systems
     if M.op=='s', nM  = feval(M.sample,nM); end
     pM = M;
     sc = 1;
    otherwise  % Use mid-point numerical differentiation
     pM = theta2m(theta+thn,Ms);
     nM = theta2m(theta-thn,Ms);
     if M.op=='s',
      % If continuous model required then sample these systems
      pM = feval(M.sample,pM);
      nM = feval(M.sample,nM);
     end
     sc = 2;
   end
   
   % Now compute numerical derivative
   dA  = (pM.ss.A-nM.ss.A)/(sc*deps);
   dC  = (pM.ss.C-nM.ss.C)/(sc*deps);
   if isB,  dB  = (pM.ss.B-nM.ss.B)/(sc*deps);   end
   if isD,  dD  = (pM.ss.D-nM.ss.D)/(sc*deps);   end
   if isK,  dK  = (pM.ss.K-nM.ss.K)/(sc*deps);   end
   if isF,  dF  = (pM.ss.F-nM.ss.F)/(sc*deps);   end
   if isG,  dG  = (pM.ss.G-nM.ss.G)/(sc*deps);   end
   if isX1, dx  = (pM.ss.X1-nM.ss.X1)/(sc*deps); end
   
  else % OK we can compute the derivative directly
   if strcmpi(M.par,'ddlc'),
    if isX1, % Need to handle the case of initial state separately
     if i<=nQp, % This gives the derivative of system matrices wrt theta
      thn = [Md.Qp(:,i);zeros(n,1)];
     else % This is derivative of initial state wrt theta
      thn = zeros(length(theta),1);
      thn(size(Md.Qp,1)+i-nQp) = 1;
     end
    else % This gives the derivative of system matrices wrt theta
     thn   = Md.Qp(:,i);
    end
    dA(:) = thn(1:n*n);         idx=n*n;
    if isB,  dB(:) = thn(idx+1:idx+n*m);   idx=idx+n*m;   end
    dC(:) = thn(idx+1:idx+p*n); idx=idx+p*n;
    if isD,  dD(:) = thn(idx+1:idx+p*m);   idx=idx+p*m;   end
    if isK,  dK(:) = thn(idx+1:idx+n*p);   idx=idx+n*p;   end
    if isF,  dF(:) = thn(idx+1:idx+n*n*m); idx=idx+n*n*m; end
    if isG,  dG(:) = thn(idx+1:idx+p*n*m); idx=idx+p*n*m; end
    if isX1, dx(:) = thn(idx+1:idx+n);     idx=idx+n;     end
   else % For all other parametrizations, we just cycle through theta
    thn = zeros(size(theta));
    thn(theta_struct(i)) = 1;
    dM  = theta2m(thn,Ms);
    dA  = dM.ss.A;
    dC  = dM.ss.C;
    if isB,  dB  = dM.ss.B;  end
    if isD,  dD  = dM.ss.D;  end
    if isK,  dK  = dM.ss.K;  end
    if isF,  dF  = dM.ss.F;  end
    if isG,  dG  = dM.ss.G;  end
    if isX1, dx  = dM.ss.X1; end
   end
  end
  
  % Run the filter to get the derivative
  if strcmpi(M.type,'ss'),
   dxh = ltitr(A-K*C,[dB-dK*D-K*dD,dK,dA-dK*C-K*dC],uyx,dx);
   dy  = xh*dC.' + dxh*Ct + u*dD.';
  else
   [dxh,dy,dpe,ukdx] = systr(0*y,bigu,A,[dB dA dF dK],C,[dD dC dG zeros(p,p)],K,F,G,dx,isD,isK,isF,isG,isX1,m);
  end
  phi(:,i) = -dy(:);
 end
 
 % Compute gradient and change phi according to cost
 switch OPT.cost,
  case 'det'
   % Modify phi to account for N
   phi=phi*scal;
   for i=1:size(phi,2),
    e        = reshape(phi(:,i),N,ny)*sqrtLam;
    phi(:,i) = e(:);
   end
   grad = 2*phi'*pe(:)/N;
   
  otherwise
   grad = 2*phi'*pe(:)/N;
 end % switch
 
 % Rotate the Jacobian and gradient into right format for argmin
 nphi = size(phi,2);
 R    = triu(qr([phi pe(:)]));
 phi  = R(1:nphi,1:nphi);
 pe   = R(1:nphi,nphi+1);
end
