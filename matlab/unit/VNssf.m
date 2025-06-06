%  Function to calculate quadratic cost
%
%  V_N(theta) = 1/N\sum_{k=1}^N[Y(w_k) - G(theta,w_k)]^2
%
%  Associated with the state-space model for observed frequency domain data:
%
%  G(theta,w_k) = C (pI - A)^{-1} B + D
%
%  The operator p can be the Z tranform variable with z=e^(j*w_k*T) (with T
%  being the sampling period in seconds), or the Laplace Transform
%  variable s with s=j*w_k.
%
%  This function is not meant to be directly called by users - instead it
%  is an auxiliary function used by foe.m for frequency domain system
%  identification.
%
%  Usage is
%
%  [cost,pe,grad,psi] = VNssf(Z,theta,OPT.n)
%
%  Z        = Frequency-Response-Function (FRF) data. Z.y(i,j,k) holds
%             the i'th output, j'th input FRF data for the k'th
%             frequency point. The frequency points are stored in Z.w.
%  theta    = [A(:);B(:);C(:);D(:)]: specification of point to calculate cost at.
%  M        = Data structure which defines the model structure which
%             is to be estimated from the data as follows:
%  OPT      = Data structure which defines options for an estimation
%             algorithm.
%
%   written by Brett Ninness,  School of EE & CS
%              Adrian Wills    University of Newcastle
%                 		           Australia.


% Copyright (C) Brett Ninness

function [cost,pe,grad,phi,map]=VNssf(Z,theta,OPT,M,div)

%Get state-space model
g=theta2m(theta,M,1);

%Find out if there is a structure that we must obey
if isfield(M,'theta_struct'),
 theta_struct = M.theta_struct;
else
 theta_struct = find(ones(length(theta),1));
end

% Compute cost associated with this state-space model
[p,m,N]=size(Z.y); n=size(g.ss.A,1);
if strcmp(M.op,'q'),
 ew=exp(j*Z.w(:)*M.T);
elseif strcmp(M.op,'s'),
 ew=j*Z.w(:);
else
 error(['OPT.op ' OPT.op ' not supported']);
end

%Get data in useful form
y=Z.y;

%Extract system matrices
A  = g.ss.A;
B  = g.ss.B;
C  = g.ss.C;
D  = g.ss.D;
K  = g.ss.K;
X1 = g.ss.X1;
Ms = M; %Save initial system

if ~isfield(OPT,'R'),
 R=ones(p*m*N,1);
else
 R=sqrt(OPT.R(:));
 OPT.cost='mse';
end
OPT.cost='mse';

%Construct matrices depending on D
isD  = ~isempty(D);
isK  = ~isempty(K);
isX1 = ~isempty(X1);

%Form (A-KC)
if isK,
 AKC = A-K*C;
else
 AKC = A;
 K   = zeros(n,p);
end
if isK && isD,
 BKD = B-K*D;
else
 BKD = B;
end
if ~isD,
 D = zeros(p,m);
end

% Get frequency response
onN  = ones(1,N); pe = zeros(N*p,m); In = eye(n,n);
if div, 
 XB   = frmimo(AKC,BKD,In,ew);
 CX   = frmimo(AKC,eye(n),C,ew);
end
pe = y - frmimo(AKC,BKD,C,ew);
if isD,
 pe = pe - D(:,:,onN);
end
if isK,
 for im=1:m,
  pe(:,im,:) = pe(:,im,:) - reshape(C*K*squeeze(y(:,im,:)),p,N);
 end
end

% Scale the prediction error according to user defined weighting
pe_saved = pe;
pe = R.*pe(:);

%Calculate cost according to selection OPT.cost
switch OPT.cost,
 case {'mse','trace'}
  cost  = real(pe'*pe)/N;
  pe    = [real(pe(:));imag(pe(:))];
 case {'det','ml'}
  pe    = reshape(pe,p,m*N);
  V     = (pe*pe')/N;
  cost  = real(log(det(V)));
  sV    = sqrtm(V)\eye(size(V));
  scal  = 1;%sqrt(det(V));
  pe    = scal*sV*pe;
  pe    = [real(pe(:));imag(pe(:))];
 otherwise
  error('Value in OPT.cost is not valid!');
end


if div,
 deps = eps^(2/3);
 
 if strcmpi(M.par,'ddlc'),
  Md  = ddlc(theta,Ms);
  nQp = size(Md.Qp,2);
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
 phi  = zeros(2*m*p*N,nth);
 dy   = zeros(p,m,N);
 dX   = zeros(n,m);
 
 %Set some default matrix sizes
 dA = zeros(n,n);
 dB = zeros(n,m);
 dC = zeros(p,n);
 dD = zeros(p,m);
 dK = zeros(n,p);
 dx = zeros(n,1);
 
 %Loop over the parameters we are allowed to vary
 for i=1:nth,
  %Compute derivative of system matrices wrt theta based on
  %parametrization (note that op='s' or M.par='grey' will result in numerical
  %derivatives being used).
  if strcmpi(M.par,'grey'),
   %Use finite difference to get derivative estimate of state-space
   %model w.r.t. theta (this will often be exact to machine eps)
   thn    = zeros(size(theta));
   thn(theta_struct(i)) = deps;
   switch lower(OPT.sysnd),
    case 'forward' %Use forward numerical differentiation
     pM = theta2m(theta+thn,Ms);
     %If continuous model required then sample these systems
     if M.op=='s', pM  = samplek(pM); end
     nM = M;
     sc = 1;
    case 'backward' %Use backward numerical differentiation
     nM = theta2m(theta-thn,Ms);
     %If continuous model required then sample these systems
     if M.op=='s', nM  = feval(M.sample,nM); end
     pM = M;
     sc = 1;
    otherwise  %Use mid-point numerical differentiation
     pM = theta2m(theta+thn,Ms);
     nM = theta2m(theta-thn,Ms);
     if M.op=='s',
      %If continuous model required then sample these systems
      pM = feval(M.sample,pM);
      nM = feval(M.sample,nM);
     end
     sc = 2;
   end
   
   %Now compute numerical derivative
   dA  = (pM.ss.A-nM.ss.A)/(sc*deps);
   dB  = (pM.ss.B-nM.ss.B)/(sc*deps);
   dC  = (pM.ss.C-nM.ss.C)/(sc*deps);
   if isD,  dD  = (pM.ss.D-nM.ss.D)/(sc*deps);   end
   if isK,  dK  = (pM.ss.K-nM.ss.K)/(sc*deps);   end
   if isX1, dx  = (pM.ss.X1-nM.ss.X1)/(sc*deps); end
   
  else %OK we can compute the derivative directly
   if strcmpi(M.par,'ddlc'),
    if isX1, %Need to handle the case of initial state separately
     if i<=nQp, %This gives the derivative of system matrices wrt theta
      thn = [Md.Qp(:,i);zeros(n,1)];
     else %This is derivative of initial state wrt theta
      thn = zeros(length(theta),1);
      thn(size(Md.Qp,1)+i-nQp) = 1;
     end
    else %This gives the derivative of system matrices wrt theta
     thn   = Md.Qp(:,i);
    end
    dA(:) = thn(1:n*n);         idx=n*n;
    dB(:) = thn(idx+1:idx+n*m); idx=idx+n*m;
    dC(:) = thn(idx+1:idx+p*n); idx=idx+p*n;
    if isD,  dD(:) = thn(idx+1:idx+p*m);   idx=idx+p*m;   end
    if isK,  dK(:) = thn(idx+1:idx+n*p);   idx=idx+n*p;   end
    if isX1, dx(:) = thn(idx+1:idx+n);     idx=idx+n;     end
   else %For all other parametrizations, we just cycle through theta
    thn = zeros(size(theta));
    thn(theta_struct(i)) = 1;
    dM  = theta2m(thn,Ms);
    dA  = dM.ss.A;
    dB  = dM.ss.B;
    dC  = dM.ss.C;
    if isD,  dD  = dM.ss.D;  end
    if isK,  dK  = dM.ss.K;  end
    if isX1, dx  = dM.ss.X1; end
   end
  end
  
  %Now run the filter to get the derivative
  dAKC = dA-dK*C-K*dC;
  dBKD = dB-dK*D-K*dD;
  dy = sysfr(N,y,XB,CX,dAKC,dBKD,dK,dC,dD,isK,isD);        
  phi(:,i) = -[real(dy(:));imag(dy(:))]; 
  if strcmpi(OPT.cost,'det'),
   e        = sV*reshape(phi(:,i),p,2*N*m);
   phi(:,i) = pe(:);
  end
 end
  
 %Now rotate the Jacobian and gradient into right format for argmin
 grad = 2*phi'*pe(:)/N;
 nphi = size(phi,2);
 R    = triu(qr([phi pe(:)]));
 phi  = R(1:nphi,1:nphi);
 pe   = R(1:nphi,nphi+1);
end
