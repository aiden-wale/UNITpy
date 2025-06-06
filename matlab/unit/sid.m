%   This function estimates a state-space model for a possibly
%   multivariable system by using a member of the `subspace
%   based' family of algorithms.
%
%   Usage is:  G = sid(Z,M,OPT);
%
%   where
%
%    Z:        Input-Output data in one of two forms.  The standard form
%              is for it to be a record with elements Z.y and Z.u, each
%              of which are matrices with number of rows equal to the
%              number of data samples, and number of columns equal (respectively)
%              to the number of outputs and the number of inputs.  On
%              the other hand, Z can be a matrix of the form Z = [y,u]
%              where it is assumed that y is a column vector of output
%              measurements and u is a matrix whose columns are the
%              input measurements; in this latter MISO models are
%              being considered.
%    M.A:      Dimension of state-space of estimated model. If set to a
%              polynomial or a matrix, only the size of this is sensed
%              for the purposes of setting the order.
%    M.delay:  Number of samples of delay to include. In the
%              case of a MIMO system, this should be a vector of delays,
%              one for each input being considered.
%    M.op:     set to 'q' for shift and 'd' for delta.  Default = 'q'.
%    M.T:      sampling period (ignored for q operator case).  Default=1
%    M.w:      vector of frequencies at which to calculate frequency
%              response of estimated model.  Specify in real frequency,
%              not normalised.  Default is 3 decades up to folding freq.
%   OPT:       Data structure which defines options for the estimation
%              algorithm as follows:
%    OPT.n:    Number of starting data points to discard to get
%              rid of initial condition effects.  Default is none.
%    OPT.bw:   Bandwidth in rad/s that model should fit over - this is
%              only used for delta operator case.  Default = 0.2/T;
%    OPT.alg:  Type of algrorithm to use.  Options are
%              OPT.alg='cca'  : Larimore's Canonical Variate Analysis.
%	       OPT.alg='n4sid': DeMoore & Overschee's N4SID (Default).
% OPT.horizon: number of lags into past that are used in forming an initial
%              basis for the state space.
%   G:         Data structure which specifies the estimated model as
%              follows:
% G.A, G.B     Matrices definining the estimated transfer function model.
% G.C, G.D     For SISO systems, these element are row vectors defining
%              co-efficients of increasing powers of M.op^-1.  For MISO,
%              they are matrices of rows, the k't row pertaining to the
%              k'th input.  For MIMO, they are 3 dim matrices with the
%              row (k,:,m) defining the transfer function from the k'th
%              input to the m'th output.
%G.ss.A,G.ss.B: [A,B,C,D] matrices/vectors defining estimated state space
%G.ss.C,G.ss.D: model.
%   G.sing:    Singular values in projection of output onto inputs.
%   G.G:       Matrix of frequency responses.  If the system has multiple
%              inputs and multpile outputs, then this matrix is 3
%              dimensional, with one `page' per output, and the i'th
%              column of the j'th page (ie G.G(:,i,j)) being the
%              frequency response from input i to ouput j.
%
%Written by Brett Ninness, School of EE & CS
%                          University of Newcastle
%                          Australia.

% Copyright (C) Brett Ninness.

function G = sid(Z,M,OPT);

% Extract inputs and outputs specified
[y,u,ny,nu,N,Z] = Z2data(Z);

% Check which parts of model structure were unspecified and set to defaults.
if ~exist('M') error('Need to specify initial model structure M!');
else
 M.type ='ss'; 
 M      = startM(Z,M);
 order  = M.nx;
end;

%Start with G = M
G = M;

% Include delays specified in model structure on inputs
for r=1:nu u(:,r) = [zeros(M.delay(r),1);u(1:N-M.delay(r),r)]; end;

% Check what options not specified explicitly by user and then set to defaults
if (nargin < 3) OPT = startOPT([]); else OPT = startOPT(OPT); end;  % Generic defaults
if ~isfield(OPT,'bw') OPT.bw =0.5/M.T;                        end;  % Defaults specific to sid.m
if ~isfield(OPT,'horizon')
 OPT.horizon=min( floor([2.5*order,(N-OPT.n+1-ny)/(2*nu+ny+2)]) );
end;
if ~isfield(OPT,'alg') OPT.alg='n4sid';
elseif strcmp(lower(OPT.alg),'sid') OPT.alg='n4sid';
elseif ~any(strcmp(lower(OPT.alg),{'n4sid','cca'})) OPT.alg='n4sid';
end;

% Forward and backward horizons are the same;
f=OPT.horizon;  p=f;

%  Form block Hankel Matrices from inputs and outputs
Y = zeros(ny*f,N-(p+f)-OPT.n+1);
U = zeros(nu*f,N-(p+f)-OPT.n+1);
Z = zeros((nu+ny)*p,N-(p+f)-OPT.n+1);

if (M.op=='q')
 for k=1:f
  Y(1+(k-1)*ny:k*ny,:) = y(OPT.n+k+p:N-f+k,:)';
  % Stack in reverse for ease of estimating state sequences from QR
  if nu>0 U(1+(f-k)*nu:((f-k)+1)*nu,:) = u(OPT.n+k+p:N-f+k,:)'; end;
 end;
 for k=1:p
  if nu>0 Z(1+(k-1)*nu:k*nu,:) = u(OPT.n+p-(k-1):N-f-k+1,:)'; end;
  Z(nu*p+1+(k-1)*ny:nu*p+k*ny,:) = y(OPT.n+p-(k-1):N-f-k+1,:)';
 end;
elseif (M.op=='d')  % Delta operator model - experimental!
 % poles = (OPT.bw/5)*exp(j*(0:1:f+p-1)*2*pi/(f+p)) - OPT.bw;
 poles = -OPT.bw*ones(1,(f+p)-1);
 J = real(poly(poles)); J = J/J(length(J));
 yf = zeros(N,ny);  uf = zeros(N,nu); zf = zeros(N,ny+nu);
 for k=1:p
  for m =1:ny
   yf(:,m) = delsimf([zeros(1,f+k-1),1,zeros(1,p-k)],J,y(:,m),M.T);
  end;
  for m =1:nu
   uf(:,m) = delsimf([zeros(1,f+k-1),1,zeros(1,p-k)],J,u(:,m),M.T);
  end;
  Z(1+(k-1)*nu:k*nu,:) = uf(OPT.n+p:N-f,:)';
  Z(nu*p+1+(k-1)*ny:nu*p+k*ny,:) = yf(OPT.n+p:N-f,:)';
 end;
 for k=1:f
  for m =1:ny
   yf(:,m) = delsimf([zeros(1,f-k),1,zeros(1,p+k-1)],J,y(:,m),M.T);
  end;
  for m =1:nu
   uf(:,m) = delsimf([zeros(1,f-k),1,zeros(1,p+k-1)],J,u(:,m),M.T);
  end;
  Y(1+(k-1)*ny:k*ny,:) = yf(OPT.n+p:N-f,:)';
  U(1+(f-k)*nu:((f-k)+1)*nu,:) = uf(OPT.n+p:N-f,:)';
 end;
end;

% Switch to ortho basis for quantities involved - much more efficient flop wise
R    = qr([U',Z',Y']); R = triu(R);
R11  = R(1:f*nu,1:f*nu);
R12  = R(1:f*nu,f*nu+1:(f+p)*nu+p*ny);
R12p = R(1:(f-1)*nu,(f-1)*nu+1:(f+p)*nu+(p+1)*ny);
R13  = R(1:f*nu,(f+p)*nu+p*ny+1:(f+p)*(nu+ny));
R22  = R(f*nu+1:(f+p)*nu+p*ny,f*nu+1:(f+p)*nu+p*ny);
R22p = R((f-1)*nu+1:(f+p)*nu+(p+1)*ny,(f-1)*nu+1:(f+p)*nu+(p+1)*ny);
R23  = R(f*nu+1:(f+p)*nu+p*ny,(f+p)*nu+p*ny+1:(f+p)*(nu+ny));
R23p = R((f-1)*nu+1:(f+p)*nu+(p+1)*ny,(f+p)*nu+(p+1)*ny+1:(f+p)*(nu+ny));
R33  = R((f+p)*nu+p*ny+1:end,(f+p)*nu+p*ny+1:(f+p)*(nu+ny));
%R33p = R((f+p)*nu+(p+1)*ny+1:end,(f+p)*nu+(p+1)*ny+1:(f+p)*(nu+ny));

% Linear regression solution wrt to orthonormal bases
beta = R23'*R22*pinv(R22'*R22);  Zp = [R12',R22'];  % Y*PI*Z'*(Z*PI*Z')^(-1) = beta
bplus= R23p'*R22p*pinv(R22p'*R22p);  Zplus=[R12p',R22p']; % For Xplus computation

% Compute weightings according to algorithm variant
if strcmp(lower(OPT.alg),'cca')
 [UU,SS,VV] = svd(R23'*R23 + R33'*R33);  % (Y_f \Pi_{U_f} Y_f^T)
 iWf = UU*sqrt(SS); Wf = pinv(iWf);      % Above to -1/2 and 1/2;
 Wp = R22'; Wplus = R22p';
elseif strcmp(lower(OPT.alg),'n4sid')
 [mm,nn] = size(beta); iWf = eye(mm,mm); Wf = iWf; Wp=Zp;
end;

% Estimate subspace spanned by columns of extended observabililty matrix O
[O,S,V] = svd(Wf*beta*Wp);  O = iWf*O(:,1:order);

% Extract estimated states - used for Q matrix estimation + par est in n4sid
Xplus = O(1:(f-1)*ny,:)\(bplus*Zplus); Xminus= O\(beta*Zp);
YY = [R13',R23',R33'];
[mm,nn] = size(YY); UU = [R11',zeros(f*nu,nn-f*nu)];
[rr,kk] = size(Xplus);  Xplus = [Xplus,zeros(rr,nn-kk)];
[rr,kk] = size(Xminus); Xminus = [Xminus,zeros(rr,nn-kk)];

% Extract estimates of A,B,C,D from the estimate of the observability
% matrix O and the estimate of the regressor beta.
if strcmp(lower(OPT.alg),'n4sid')
 % If using n4sid method, then we are regressing on estimated states
 %TH = [Xplus;YY(1:ny,:)]/[Xminus;UU((f-1)*nu+1 :f*nu,:)];
 TH = [Xplus;YY(1:ny,:)]*pinv([Xminus;UU((f-1)*nu+1 :f*nu,:)]);
 G.ss.A = TH(1:order,1:order);
 G.ss.B = TH(1:order,order+1:order+nu);
 G.ss.C = TH(order+1:order+ny,1:order);
 G.ss.D = TH(order+1:order+ny,order+1:order+nu);
elseif strcmp(lower(OPT.alg),'cca')
 %  estimate A and C using Kung's algorithm.
 [m,n] = size(O);
 C = O(1:ny,:);  A = O(1:m-ny,:)\O(ny+1:m,:);
 G.ss.A = A; G.ss.C = C;
 
 %  Estimate B and D using least squares
 if nu>0  % Only do this if there is an exogenous input
  yy = reshape(y',N*ny,1);
  for k =1:nu
   phid(:,(k-1)*ny+1:k*ny)=kron(u(:,k),eye(ny,ny));
   for m = 1:order
    [num,den] = ss2tf(A,[zeros(m-1,1);1;zeros(order-m,1)],C,zeros(ny,1));
    for n =1:ny
     indexs = kron(ones(N,1),[zeros(1,n-1),1,zeros(1,ny-n)]);
     if (M.op=='q')
      zz = filter(num(n,:),den,u(:,k));
     else
      zz = delsimf(num(n,:),den,u(:,k),M.T);
     end;
     phib(logical(indexs),(k-1)*order+m) = zz;
    end;
   end;
  end;
  PHI = [phid,phib]; th = PHI\yy;
  Dvec = th(1:nu*ny); Bvec = th(nu*ny+1:length(th));
  G.ss.D = reshape(Dvec,ny,nu); G.ss.B = reshape(Bvec,order,nu);
 else
  G.ss.B=[]; G.ss.D=[];
 end;
 
 % Put in stacked form for subsequent Q and R matrix estimation
 TH(1:order,1:order) = G.ss.A; TH(1:order,order+1:order+nu) = G.ss.B;
 TH(order+1:order+ny,1:order) = G.ss.C; TH(order+1:order+ny,order+1:order+nu) = G.ss.D;
end;

% Compute state and noise covariances from residuals
ew = [Xplus;YY(1:ny,:)] - TH*[Xminus;UU((f-1)*nu+1 :f*nu,:)];
w=ew(1:order,:); e = ew(order+1:order+ny,:); nn=size(w,2);
G.ss.Q = w*w'/nn;  G.ss.R = e*e'/nn; G.ss.S = w*e'/nn;

% Convert from state space to transfer function form
% firstly remove G.A and G.B in case they came with the initial model
if isfield(G,'A'), G = rmfield(G,'A'); end
if isfield(G,'B'), G = rmfield(G,'B'); end
for k=1:nu
 for m=1:ny
  [G.B(k,:,m),G.A(k,:,m)] = ss2tf(G.ss.A,G.ss.B,G.ss.C(m,:),G.ss.D(m,:),k);
 end;
end;

% Return final estimate in innovations form - suppress warnings/errors from dare
try
 [P,dum,dum,rep] = dare(G.ss.A',G.ss.C',G.ss.Q,G.ss.R,G.ss.S);
 V=G.ss.C*P*G.ss.C'+G.ss.R; G.ss.K=(G.ss.A*P*G.ss.C'+G.ss.S)*inv(V)';
catch
 V=eye(size(G.ss.A)); G.ss.K = zeros(order,ny);
end;

% Convert from state space to transfer function form
G = sstotf(G);

% Fill in last remaining details about the estimated model structure.
G.T = M.T; 
G.op = M.op; 
G.w = M.w;
G.sing = diag(S(1:f,1:f));  
G.type = 'ss';  
G.delay = M.delay;
if isfield(M,'in'),
 G.in = M.in;
end
if isfield(M,'out'),
 G.out = M.out;
end

% Fill in bilinear info
G.ss.F  = [];
G.ss.G  = [];
G.ss.X1 = zeros(order,1);

% Add legend for prospective plotting
if strcmp(lower(OPT.alg),'n4sid')
 G.disp.legend=['Estimated SS model via N4SID'];
elseif strcmp(lower(OPT.alg),'cca')
 xsG.disp.legend=['Estimated SS model via CCA'];
end;

% Add in estimation of innovations variance
G.var = trace(V);

%G.alg='sid'; % Record that subspace method was used
G.alg=OPT.alg; % Record that subspace method was used
G.OPT = OPT;

%Record that VNss should be used to compute prediction erros for validation
G.costfcn = 'VNss';

%Finally, make sure that M.theta reflects the model
G.theta = [G.ss.A(:);G.ss.B(:);G.ss.C(:)];
if G.estD,  G.theta = [G.theta;G.ss.D(:)]; end
if G.estK,  G.theta = [G.theta;G.ss.K(:)]; end
if G.estX1, G.theta = [G.theta;G.ss.X1];   end

