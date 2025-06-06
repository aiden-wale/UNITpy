% Square-root implementation of Kalman Filter/Smoother as per
%
% Gibson and Ninness, "RobustMaximum-Likelihood Estimation of Multivariable
% Dynamic Systems", Automatica, , 41(10):1667?1682, October 2005.
%
% Model is assumed to be of the form,
%
% x(t+1) = A(t)x(t) + B(t)u(t) + w(t),   [w(t)]    (    [ Q(t)    S(t) ] )
%                                        [    ] ~ N( 0, [              ] )
%   y(t) = C(t)x(t) + D(t)u(t) + v(t),   [v(t)]    (    [ S^T(t)  R(t) ] )
%
%  Usage is
%
%      G = rksqrtv(Z,M,OPT);
%
% Where
%
%             Z: the data structure containing the measured outputs in Z.y and 
%                possibly the measured inputs in Z.u. If the number of data
%                points is N, the number of outputs is p, and the number of
%                inputs is m, then Z.y is an N x p matrix and Z.u is an N x m
%                matrix.
%
%  M.ss.A,B,C,D: 3D time-varying system matrices; for each matrix, it is
%                assumed that the time index is the third dimension, e.g.
%                A(:,:,t) is the state transition matrix at time t
%
%    M.ss.Q,S,R: 3D time-varying covariance matrices for process and measurement 
%                noise, respectively. Again, the third dimension is time
%                dependent.
%
%    M.ss.X1,P1: Initial state mean (X1) and its covariance matrix (P1),
%                respectively.
%
%             G: returned structure with the following fields
%
%          
%          G.xp: predicted states, i.e. E[x(t) | y_1,..,y_{t-1}]
%          G.xf: filtered states, i.e. E[x(t) | y_1,..,y_{t}]
%          G.xs: smoothed states, i.e. E[x(t) | y_1,..,y_{N}]
%
%          G.Pp: squareroot of predicted state covariance matrix
%                  Pp(:,:,t) = E{x(t)*x(t)' |  y_1,..,y_{t-1}}
%          G.Pf: squareroot of filtered state covariance matrix
%                  Pf(:,:,t) = E{x(t)*x(t)' |  y_1,..,y_{t}}
%          G.Ps: squareroot of smoothed state covariance matrix, i.e.
%                  Ps(:,:,t) = E{x(t)*x(t)' |  y_1,..,y_{N}}
%          G.Ms: cross covariance between x(t+1) and x(t), i.e.
%                  Ms(:,:,t) = E{x(t+1)*x(t)' | y_1,..,y_{N}}
%          
%          G.yp: predicted output estimate
%          G.yf: filtered output estimate
%          G.ys: smoothed output estimate
%          G.pe: prediction error
%          G.fe: filter error
%          G.se: smoother error
%
%          G.LL: negative log-likelihood
%
%
%           OPT: a structure containing algorithm options. 
% OPT.smoothing: if equal to 1, then smoothing will be performed in
%                addition to filtering.
%
%
%   written by Brett Ninness: School of EE & CS
%              Adrian Wills   University of Newcastle
%        	                  Australia.
%

% Copyright Brett Ninness

function G = rksqrtv(Z,M,OPT)

% Extract data and sizes
y = Z.y;
u = Z.u;
p = size(y,1);
m = size(u,1);
N = length(y);

%If we have no measurements then we cant filter or smooth
if p<1,
 error('According to the data, there are no outputs. Nothing to do.')
end

%Check to see if we are smoothing as well as filtering
if ~isfield(OPT,'smoothing'),
 OPT.smoothing = 0;  %Default is filtering only
end

%Check to see if we have been called by the EM function
if ~isfield(OPT,'ksem'),
 OPT.ksem = 0;  %Default is NO EM H matrix
end

%Need to determine what type of state space system we are dealing with and
%if all the necessary matrices are present
if strfind(M.type,'bilin'),
 bilin = 1;
else
 bilin = 0;
end

% Let G = M
G = M;

%Extract system matrices and noise model
A  = M.ss.A;
n  = size(A,1);
B  = M.ss.B;
C  = M.ss.C;
D  = M.ss.D;
Q  = M.ss.Q;
R  = M.ss.R;
S  = M.ss.S;
P1 = M.ss.P1;
X1 = M.ss.X1;

if strfind(M.type,'bilin'),
 if ~isfield(M.ss,'F'), M.ss.F = zeros(n,n*max(1,m)); end
 if ~isfield(M.ss,'G'), M.ss.G = zeros(p,n*max(1,m)); end
 
 % If the system is a bilinear system then adjust the A and C matrices
 A = zeros(n,n,N);
 C = zeros(p,n,N);
 for t=1:N,
  uu       = kronaw(u(:,t),eye(max(1,m)));
  A(:,:,t) = M.ss.A + M.ss.F*uu;
  C(:,:,t) = M.ss.C + M.ss.G*uu;
 end
end


%The matrix H will store the expected value of
%
%   [      x(t)     ]   [      x(t)     ]'
%   [      u(t)     ]   [      u(t)     ] 
%   [kron(u(t),x(t))] * [kron(u(t),x(t))]
%   [     x(t+1)    ]   [     x(t+1)    ] 
%   [      y(t)     ]   [      y(t)     ]
%
% This will be returned from the Kalman smoother routine if OPT.ksem == 1
if OPT.ksem,
 H    = zeros(n+m+bilin*n*m+n+p);
 ix   = 1 : n;
 iu   = n+1 : n+m;
 iukx = n+m+1 : n+m+bilin*n*m;
 ix1  = n+m+bilin*n*m+1 : n+m+bilin*n*m+n;
 iy   = n+m+bilin*n*m+n+1 : n+m+bilin*n*m+n+p;
end

%Make some room for mean and covariance
xp = zeros(n,N+1);
xf = zeros(n,N);

Pp = zeros(n,n,N+1);
Pf = zeros(n,n,N);
Ri = zeros(p,p,N);
K  = zeros(n,p,N);

yp = zeros(p,N);
yf = zeros(p,N);
pe = zeros(p,N);
fe = zeros(p,N);

%Initialise the predicted mean and cov
xp(:,1)   = X1;
Pp(:,:,1) = rchol(P1);

%Make some room
LL  = 0;
PE  = 0;
R1  = zeros(n+p);
R2  = zeros(n+p);
R3  = zeros(2*n,n);
R4  = zeros(2*n,n);

%Are the system matrices (indicated by last letter) time-varying
tmp = size(A,3);
if tmp > 1 && tmp ~= N,
 error(['The A matrix is specified as a time varying matrix but it does not have the correct number of time entries, i.e. 3rd dimension is ' num2str(tmp) ' but should be ' num2str(N) '.']);
end
tva = tmp>1;
tmp = size(B,3);
if tmp > 1 && tmp ~= N,
 error(['The B matrix is specified as a time varying matrix but it does not have the correct number of time entries, i.e. 3rd dimension is ' num2str(tmp) ' but should be ' num2str(N) '.']);
end
tvb = tmp>1;
tmp = size(C,3);
if tmp > 1 && tmp ~= N,
 error(['The C matrix is specified as a time varying matrix but it does not have the correct number of time entries, i.e. 3rd dimension is ' num2str(tmp) ' but should be ' num2str(N) '.']);
end
tvc = tmp>1;
tmp = size(D,3);
if tmp > 1 && tmp ~= N,
 error(['The D matrix is specified as a time varying matrix but it does not have the correct number of time entries, i.e. 3rd dimension is ' num2str(tmp) ' but should be ' num2str(N) '.']);
end
tvd = tmp>1;
tmp = size(Q,3);
if tmp > 1 && tmp ~= N,
 error(['The Q matrix is specified as a time varying matrix but it does not have the correct number of time entries, i.e. 3rd dimension is ' num2str(tmp) ' but should be ' num2str(N) '.']);
end
tvq = tmp>1;
tmp = size(S,3);
if tmp > 1 && tmp ~= N,
 error(['The S matrix is specified as a time varying matrix but it does not have the correct number of time entries, i.e. 3rd dimension is ' num2str(tmp) ' but should be ' num2str(N) '.']);
end
tvs = tmp>1;
tmp = size(R,3);
if tmp > 1 && tmp ~= N,
 error(['The R matrix is specified as a time varying matrix but it does not have the correct number of time entries, i.e. 3rd dimension is ' num2str(tmp) ' but should be ' num2str(N) '.']);
end
tvr = tmp>1;

%If S is nonzero and (Q or S or R) are time varying, then Q and A and B will also 
% be time varying (if they are not already), so we need to make some extra room to store them
normS = norm(S(:));
if any([tvq tvs tvr]) && normS>0.0,
 if ~tva,
  A   = A(:,:,ones(1,N));
  tva = 1;
 end
 if ~tbv,
  B   = B(:,:,ones(1,N));
  tvb = 1;
 end
 if ~tvq,
  Q   = Q(:,:,ones(1,N));
  tvq = 1;
 end
end

%Current index of time for time-varying systems (will remain at 1 for
%time-invariant systems)
tta = 1;
ttb = 1;
ttc = 1;
ttd = 1;
ttq = 1;
tts = 1;
ttr = 1;

%Run forward Kalman filter loop
for t=1:N,
 % Robust computation of Q-S*inv(R)*S' that preserves symmetry and
 % non-negativity
 if any([tvq tvs tvr]) || t==1,
  X                   = triu(rchol([R(:,:,ttr) S(:,:,tts)'; S(:,:,tts) Q(:,:,ttq)]));
  R(:,:,ttr)          = X(1:p,1:p);
  Q(:,:,ttq)          = X(p+1:p+n,p+1:p+n);
  if normS>0.0,
   SR1             = (X(1:p,p+1:p+n)')/(X(1:p,1:p)');
  else
   SR1             = zeros(n,p);
  end
 end
 if normS>0.0,
  if any([tva tvc tvq tvs tvr]) || t==1,
   A(:,:,tta)      = A(:,:,tta)-SR1*C(:,:,ttc);
  end
  if any([tvb tvd tvq tvs tvr]) || t==1,
   B(:,:,ttb)      = B(:,:,ttb)-SR1*D(:,:,ttd);
  end
 end
 R1(1:p,1:p)         = X(1:p,1:p);
 R1(p+1:end,1:p)     = Pp(:,:,t)*C(:,:,ttc)';
 R1(p+1:end,p+1:end) = Pp(:,:,t);
 R2                  = triu(qr(R1)); 
 Ri(:,:,t)           = R2(1:p,1:p);
 K(:,:,t)            = R2(1:p,p+1:p+n)'/(R2(1:p,1:p)');
 Pf(:,:,t)           = R2(p+1:p+n,p+1:p+n);
 yp(:,t)             = C(:,:,ttc)*xp(:,t) + D(:,:,ttd)*u(:,t);
 pe(:,t)             = y(:,t) - yp(:,t); 
 PE                  = PE + pe(:,t)'*pe(:,t);
 %Riep                = rfbs(pe(:,t),R2(1:p,1:p)',1);
 Riep                = R2(1:p,1:p)'\pe(:,t);
 LL                  = LL + Riep(:)'*Riep(:) + 2*sum(log(abs(diag(R2(1:p,1:p)))));
 xf(:,t)             = xp(:,t) + R2(1:p,p+1:end)'*Riep; 
 yf(:,t)             = C(:,:,ttc)*xf(:,t) + D(:,:,ttd)*u(:,t);
 fe(:,t)             = y(:,t) - yf(:,t); 
 xp(:,t+1)           = A(:,:,tta)*xf(:,t) + B(:,:,ttb)*u(:,t)  + SR1*y(:,t); 
 R3(1:n,:)           = Pf(:,:,t)*A(:,:,tta)';
 R3(n+1:end,:)       = X(p+1:p+n,p+1:p+n); 
 R4                  = triu(qr(R3,0));
 Pp(:,:,t+1)         = R4(1:n,:);
 
 %Update time-varying time index
 if t<N,
  tta = tta + tva*1;
  ttb = ttb + tvb*1;
  ttc = ttc + tvc*1;
  ttd = ttd + tvd*1;
  ttq = ttq + tvq*1;
  tts = tts + tvs*1;
  ttr = ttr + tvr*1;
 end
end

%Save parts that belong to filtering
G.ss.xf  = xf;
G.ss.Pf  = Pf;
G.ss.xp  = xp;
G.ss.Pp  = Pp;
G.ss.K   = K;
G.ss.Ri  = Ri;
G.yf     = yf;
G.yp     = yp;
G.LL     = LL;
G.PE     = PE;
G.pe     = pe;
G.fe     = fe;

if OPT.smoothing,
 %Make some room
 xs = zeros(n,N+1);
 Ps = zeros(n,n,N+1);
 Ms = zeros(n,n,N);
 ys = zeros(p,N);
 se = zeros(p,N);
 R5 = zeros(3*n,2*n);
 R6 = zeros(3*n,2*n);
 
 %Now run backward filtering (smoother stage)
 xs(:,N+1)   = xp(:,N+1);
 xs(:,N)     = xf(:,N);
 Ps(:,:,N+1) = Pp(:,:,N+1);
 Ps(:,:,N)   = Pf(:,:,N);
 APf         = A(:,:,tta)*Pf(:,:,N)';
 AP          = APf*Pf(:,:,N);
 Ms(:,:,N)   = AP;
 KCN         = eye(n)-K(:,:,N)*C(:,:,ttc);
 
 %Update the H matrix if necessary
 if OPT.ksem,
  if bilin,
   ukx    = kron(u(:,N),xs(:,N));
   v      = [xs(:,N);u(:,N);ukx;xs(:,N+1);y(:,N)];
  else
   v      = [xs(:,N);u(:,N);xs(:,N+1);y(:,N)];
  end
  Pt         = Ps(:,:,N)'*Ps(:,:,N);
  H          = v*v';
  H(ix,ix)   = H(ix,ix)   + Pt;
  H(ix,ix1)  = H(ix,ix1)  + Ms(:,:,N)';
  H(ix1,ix1) = H(ix1,ix1) + Ps(:,:,N+1)'*Ps(:,:,N+1);
  if bilin,
   H(ix,iukx)   = H(ix,iukx)   + kron(u(:,N)',Pt);
   H(iukx,iukx) = H(iukx,iukx) + kron(u(:,N)*u(:,N)',Pt);
   H(iukx,ix1)  = H(iukx,ix1)  + kron(u(:,N),Ms(:,:,N)');
  end
 end
 
 for t=N-1:-1:1,
  %Adjust time-varying time index backwards
  tta = tta - tva*1;
  ttb = ttb - tvb*1;
  ttc = ttc - tvc*1;
  ttd = ttd - tvd*1;
  ttq = ttq - tvq*1;
  tts = tts - tvs*1;
  ttr = ttr - tvr*1;
  APf                   = A(:,:,tta)*Pf(:,:,t)';
  AP                    = APf*Pf(:,:,t);
  Jt                    = (AP'/Pp(:,:,t+1))/Pp(:,:,t+1)';
  R5(1:n,1:n)           = APf';
  R5(1:n,n+1:end)       = Pf(:,:,t);
  R5(n+1:2*n,1:n)       = Q(:,:,ttq);
  R5(2*n+1:end,n+1:end) = Ps(:,:,t+1)*Jt';
  R6                    = triu(qr(R5,0));
  Ps(:,:,t)             = R6(n+1:2*n,n+1:2*n);
  xs(:,t)               = xf(:,t) + Jt*(xs(:,t+1)-xp(:,t+1));
  ys(:,t)               = C(:,:,ttc)*xs(:,t) + D(:,:,ttd)*u(:,t);
  se(:,t)               = y(:,t) - ys(:,t);
  if t==N-1,
   Ms(:,:,t)         = KCN*AP;
  else
   Ms(:,:,t)         = (Pf(:,:,t+1)'*Pf(:,:,t+1) + Jtp1*APp1)*Jt';
  end
  Jtp1                  = Jt;
  APp1                  = Ms(:,:,t) - AP;
  
  %Update the H matrix if necessary
  if OPT.ksem,
   if bilin,
    ukx    = kron(u(:,t),xs(:,t));
    v      = [xs(:,t);u(:,t);ukx;xs(:,t+1);y(:,t)];
   else
    v      = [xs(:,t);u(:,t);xs(:,t+1);y(:,t)];
   end
   Pt         = Ps(:,:,t)'*Ps(:,:,t);
   H          = H          + v*v';
   H(ix,ix)   = H(ix,ix)   + Pt;
   H(ix,ix1)  = H(ix,ix1)  + Ms(:,:,t)';
   H(ix1,ix1) = H(ix1,ix1) + Ps(:,:,t+1)'*Ps(:,:,t+1);
   if bilin,
    H(ix,iukx)   = H(ix,iukx)   + kron(u(:,t)',Pt);
    H(iukx,iukx) = H(iukx,iukx) + kron(u(:,t)*u(:,t)',Pt);
    H(iukx,ix1)  = H(iukx,ix1)  + kron(u(:,t),Ms(:,:,t)');
   end
  end
 end
 
 %Now record parts that belong to smoothing
 G.ss.xs  = xs;
 G.ss.Ps  = Ps;
 G.ss.Ms  = Ms;
 G.ys     = ys;
 G.se     = se;
end

%If necessary we need to make H symmetric
if OPT.ksem,
 idxu = find(triu(ones(size(H))));
 G.H  = H'; 
 G.H(idxu) = H(idxu);
end

return;

%--------------------------------------------------------------------------
%
%  AUXILIARY FUNCTIONS
%
%--------------------------------------------------------------------------
% Function to compute Cholesky factor robustly
function [A] = rchol(A)
A = triu(A); n = size(A,1); tol = n*eps;
if A(1,1) <= tol,
 A(1,1:n) = 0;
else
 A(1,1:n) = A(1,1:n)/sqrt(A(1,1));
end
for j=2:n,
 A(j,j:n) = A(j,j:n) - A(1:j-1,j)'*A(1:j-1,j:n);
 if A(j,j) <= tol,
  A(j,j:n) = 0;
 else
  A(j,j:n) = A(j,j:n)/sqrt(A(j,j));
 end
end


% Function that performs robust forward or backward substitution
function X = rfbs(B,A,uplow)
%uplow = 1 for lower triangular, 0 for upper
deps=100*eps; [n,m]=size(B); X=zeros(n,m);
if uplow
 for i=1:m,
  if abs(A(1,1))>deps,
   X(1,i)=B(1,i)/A(1,1);
  end
  for j=2:n,
   if abs(A(j,j))>deps,
    X(j,i)=(B(j,i)-A(j,1:j-1)*X(1:j-1,i))/A(j,j);
   end
  end
 end
else
 for i=1:m,
  if abs(A(n,n))>deps,
   X(n,i)=B(n,i)/A(n,n);
  end
  for j=n-1:-1:1,
   if abs(A(j,j))>deps,
    X(j,i)=(B(j,i)-A(j,j+1:n)*X(j+1:n,i))/A(j,j);
   end
  end
 end
end