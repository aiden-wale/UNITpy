%  Computes ARX model:
%
%  A(p)F(w_s) = B(p)cos(2*pi*w_s) + V(w_s)
%
%  from frequency domain data observations F(w_s) that are corrupted by
%  additive noise V(w_s).  The operator p can be the Z tranform variable
%  with z=e^(j*w_s*T) (with T being the sampling period in seconds), the
%  Euler differencing (delta) operator d = (q-1)/T with d =
%  (e^(j*w_s*T)-1)/T or the Laplace Transform variable s with s=j*w.  A
%  quadratic (least squares) loss criterion is used.
%
%
%   Usage is:  G = farx(Z,M,OPT);
%
%   where
%
%  Z          = observed frequency response data [F(:),w(:)] where
%               plot(w,abs(F)) should plot the measured  frequency
%               response.   Units for w are real *not* normalised freq.
%  M          = Data structure which defines the model structure which
%               is to be estimated from the data as follows:
%   M.A       = Number of poles to be estimated in denominator - which is
%               then set as equal to # of zeros to be estimated in numerator.
%   M.op      = set to 'q' for shift, 'd' for delta, 's' for Laplace
%               Default = 's'.
%   M.T       = sampling period in s. (Ignored for q case) Default = 1;
%  OPT        = Data structure which defines options for the estimation
%               algorithm as follows:
%   OPT.basis = only applicable for 's' operator models, and selects either
%               Chebychev ('cheby') or Laguerre ('ortho') orthonormal bases,
%               or normal non-orthonormal polynomial ('polyb') basis.
%               Default is 'ortho'.
%   OPT.W     = Vector of same dimension as w that specifies a
%               frequency weighting for the least squares fit.  That is,
%               plot(w,W) should give a graphical interpretation of the
%               weighting.  The default is a flat (unprejudiced) weighting.
%  G          = Data structure which specifies the estimated model as
%               follows:
%   G.B/G.A   = estimated transfer function for model of dynamics.
%   G.G       = Frequency response of estimated model for dynamics.
%   G.th      = Estimated Parameter vector from which G.B, G.A are formed.
%
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                     		  Australia.

% Copyright (C) Brett Ninness.

function G = farx(Z,M,OPT);

% Pass input data to output, then ovewrite as necessary below
G = M;

% Extract out relevant vectors from input data
[F,w,ny,nu,Ny] = Z2data(Z); F=squeeze(F); F=F(:); wmax = 1;%max(w);

% Make input equal to the identity for each w(k) if not supplied by user
if ~isfield(Z,'u'),
 Z.u = zeros(1,1,Ny);
 for k=1:Ny, Z.u(:,:,k) = 1; end
elseif isempty(Z.u),
 Z.u = zeros(1,1,Ny);
 for k=1:Ny, Z.u(:,:,k) = 1; end
end
U = squeeze(Z.u);
U = U(:);

% Check what options not specified explicitly by user and then set to
% defaults
if ~exist('OPT')
 OPT.basis = 'ortho';  OPT.W = ones(size(F));
else
 if ~isfield(OPT,'basis') OPT.basis = 'ortho';  end;
 if ~isfield(OPT,'W') OPT.W=ones(size(F));      end;
 if (length(OPT.basis) ~=5)
  error('Not a recognised basis from: ortho,cheby,polyb'); end;
end;

% Check which parts of model structure were unspecified and set to defaults.

if ~exist('M') error('Need to specify initial model structure M!');
else
 if ~isfield(M,'op')    M.op='s';   end;
 if ~isfield(M,'T')     M.T=1;      end;
 if ~isfield(M,'B')     M.B=M.A;    end;
 if ~isfield(M,'delay') M.delay=0;  end;
 if ~isfield(M,'w')     M.w= logspace(log10(pi/M.T/1000),log10(pi/M.T)); end;
 M.A = M.A(:); M.B = M.B(:);
end;

% Figure out the numerator and denominator orders
n  = M.nA+1;
m  = M.nB+1;
mn = max(m,n);

normw = 0;  %  Default is no frequency normalisation

%  Decide on frequency domain version of operator specified.

if (M.op=='q')
 z = exp(j*w*M.T);
elseif (M.op == 's')
 if (OPT.basis=='ortho')
  z = j*w;
 else
  z = j*w/wmax;
  normw = 0;
 end
elseif (M.op == 'd')
 z =  (exp(j*w*M.T)-ones(size(w)))/M.T;
end

%  Generate regressors dependent on operator and basis
zz  = zeros(length(w),mn);
PHI = zeros(length(w),m+n-1);
X   = eye(mn,mn);
%  Then use it to generate regressors
for k=1:mn
 zz(:,k) = polyval([X(mn-k+1,:)],z);
end
zz = fliplr(zz);  % order L->R from lowest order to highest order poly

%  First get regressors on u.
PHI(:,1:m)     =  (U(:)*ones(1,m)).*zz(:,end-m+1:end);
%  Then get regressors on y.
PHI(:,m+1:end) = -(F(:)*ones(1,n-1)).*zz(:,end-n+2:end); 

%  Incorporate any frequency domain weighting that was specified
F   = (OPT.W(:)).*F;
PHI = PHI.*(sqrt(OPT.W(:))*ones(1,m+n-1));

%  Can now calculate least squares estimate
FF = F.*zz(:,1); 
ph = [real(PHI);imag(PHI)]; 
ff = [real(FF);imag(FF)];
th = ph\ff; 
pe = ff - ph*th;

%  Extract A and B polynomials from theta parameter vector
%  Revert basis to normal polynomial one.
G.B = th(1:m)';
G.A = [1 th(m+1:end)'];

A1 = G.A(1); G.B = G.B/G.A(1); G.A = G.A/G.A(1); % Will need original G.A(1) later.

% Stack estimation results into output structure.
G.T = M.T; G.w = M.w; G.op = M.op; G.th = th;
G.X = X; G.delay=M.delay; G.C=[]; G.D=[];
G.PHI = PHI; G.type ='farx'; G.in.type = 'linear'; G.out.type = 'linear';
G.nB = M.nB; G.nA = M.nA;
G.nu = M.nu;
G.ny = M.ny;

% Add legend for prospective plotting
G.disp.legend=['Estimated ',G.type,' model'];

% Get Estimate of white noise variance by sample variance of residuals
G.var = pe'*pe/length(pe);

% Use it to generate estimate of parameter covariance matrix
G.P = G.var*pinv(ph'*ph);

% That variance is with respect to a particular basis; undo this
PB = G.P(1:m,1:m);  %  Easy for B
% For A, have to take into account that one element is fixed, and hence
% zero variance
PA  = [zeros(1,n);zeros(n-1,1),G.P(m+1:end,m+1:end)];
PB  = PB/(A1^2); 
PA  = PA/(A1^2); 
PA  = PA(2:n,2:n);

% Use above variances to Work out standard deviations on estimated parameters
G.SD.B = real(sqrt(abs(diag(PB)))); 
G.SD.A = real(sqrt(abs(diag(PA))));
G.SD.B = G.SD.B(:)'; 
G.SD.A = G.SD.A(:)';

G.alg='block'; % Record that block solution was used
