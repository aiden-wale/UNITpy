%  This routine takes a record of input-output data [y,u] and fits
%  an ARX (equation error) model of the form
%
%  A(p)y_t = B(p)u_{t-delay} + e_t
%
%  to [u,y] where p may be either the forward shift operator q or the Euler
%  differencing (delta) operator d = (q-1)/T where T is the sampling
%  period. The polynomials have the following structure in p^{-1}
%
%  A(p) = 1.0 + a_1 p^{-1} + a_2 p^{-2} + ... + a_n p^{-n}
%  B(p) = b_0 + b_1 p^{-1} + b_2 p^{-2} + ... + b_m p^{-m}
%
%  Usage is:
%
%  G = barx(Z,M,OPT);
%
%  where
%
%   Z:         Input-Output data in one of two forms.  The standard form
%              is for it to be a record with elements Z.y and Z.u, each
%              of which are matrices with number of rows equal to the
%              number of data samples, and number of columns equal (respectively)
%              to the number of outputs and the number of inputs.  On
%              the other hand, Z can be a matrix of the form Z = [y,u]
%              where it is assumed that y is a column vector of output
%              measurements and u is a matrix whose columns are the
%              input measurements; in this latter MISO models are
%              being considered.
%
%              In the special case that via the above only y is
%              specified, then it is assumed that an auto-regressive (AR)
%              model of the form A(p)y_t = e_t is required.
%
% M          = Data structure which defines the model structure which
%              is to be estimated from the data as follows:
%   M.A,M.B  = Number of poles $n$ and zeros $m$ to be estimated in
%              numerator and denominator poly's A(p) and B(p).
%   M.delay  = Number of samples of delay to include (see above model).  In the
%              case of a MISO system, this should be a vector of delays,
%              one for each input being considered.
%   M.op     = set to 'q' for shift and 'd' for delta.  Default = 'q'.
%   M.T      = sampling period in seconds. (Ignored for q case) Default = 1;
%   M.J      = So-called `observer polynomial' associated with delta operator
%              implementation of ARX estimation.  Ignored for op='q'.
%              to be normalised from the left.  Default = 1;
% OPT        = Data structure which defines options for the estimation
%              algorithm as follows:
%   OPT.n    = number of starting data points to discard to get
%              rid of initial condition effects.  Default is none.
%   OPT.fast = When set to 1, then makes algorithm run maximally fast by
%              avoiding the calculation of error bounds.  Default is 0.
%   OPT.filt = When set to 1, then algorithm only does filtering to
%              generate G.phi - no least squares estimation is done.
%   OPT.alg    This structure sets the algorithm used to find the estimate
%	       and also sets any parameters that are associated with
%	       the algorithm. The components of this structure are as follows:
%   OPT.alg.type - determines the algorithm type used.  It may be set as:
%	      'block' - the solution for the least squares
%			estimate is found in block form (default).
%	      'rls'   - the solution is found recursively
%			the via recursive least squares algorithm.
%	      'ct'    - the solution is found recursively via recursive
%			least squares with a contant trace forced
%			on the covariance matrix.
%	      'lms'   - the solution is found recursively via
%			least-mean-square (gradient descent) algorithm.
%	      'kf'    - the solution is found recursively via
%			a Kalman Filtering algorithm.
%   OPT.alg.P      - Initialisation of covariance matrix, default is P=10*I
%   OPT.alg.mu     - LMS gain, default is mu=0.001;
%   OPT.alg.lambda - RLS alg. `forgetting factor'. Default is lambda = 1;
%   OPT.alg.R      - Measurement noise variance used by kf alg. Default = 1;
%   OPT.alg.Q      - Parameter variance used by kf alg. Default = 0.1*I;
%   OPT.alg.th     - Initial parameter vector value.  Default = 0.
%
% G          = Data structure which specifies the estimated model as
%              follows:
%   G.B/G.A  = estimated transfer function for model of dynamics.
%   G.C/G.C  = estimated transfer function for model of noise.
%   G.G      = Frequency response of estimated model for dynamics.
%   G.H      = Frequency response of estimated model for noise colouring.
%   G.Ge     = Matrix specifying 95% confidence regions for estimated
%              frequency response Ghat.  They may be plotted by using either
%              of the commands `shownyq(G)' or `showbode(G)'.
%   G.P      = Covariance Matrix for Estimated Parameters.
%   G.th     = Parameter estimates as a column vector.
%   G.phi    = Regression matrix used to compute G.th
%
%   If a recursive solution algorithm is used then also available is:
%
%   G.th_hist History of how parameter space estimates evolved.
%   G.pe      Time history of how prediction error evolved.
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%        		              Australia.

% Copyright (C) Brett Ninness, Adrian Wills

function G = barx(Z,M,OPT)

% Extract input and output from data matrix
[y,u,ny,nu,Ny] = Z2data(Z);

% Unspecified parts of OPT -> defaults
if ~exist('OPT','var'),
 OPT = startOPT([]);
else
 OPT = startOPT(OPT);
end
if (OPT.n>=Ny),
 error('Cannot OPT.n larger than height of Z!');
end
if ~isfield(OPT,'alg'),
 OPT.alg.type = 'block';
end
if ~isfield(OPT.alg,'type'),
 OPT = rmfield(OPT,'alg');
 OPT.alg.type = 'block';
end

% Unspecified parts of M -> defaults
if ~exist('M','var'),
 error('Need to specify initial model structure M!');
else
 M = startM(Z,M);
end

%Extract number of inputs and outputs from model
nu = M.nu;
ny = M.ny;


% Include delays specified in model structure on inputs
for r=1:nu,
 u(:,r) = [zeros(M.delay(r),1);u(1:Ny-M.delay(r),r)];
end

% Compute number of numerator (and possibly denominator) parameters to be estimated
if nu>0,  % Not the AR case
 m = sum(M.nB)+nu;
else      % The AR case
 m = 0;   % Number of B(q) terms to estimate is zero
end

% For MISO ARX models, all denominator orders must be the same - set as
% the largest one specified by the user.
n    = floor(max(M.nA));
M.nA = ones(max(ny,nu),1)*n;  % max is sneaky way to handle AR case

%  Set up observer poly - only used for delta operator case
if (M.op=='d'),
 Jord = max(n,max(M.nB));
 if ~isfield(M,'J'),  % Choose default poly so delta gives identical regressors to shift
  M.J = poly((-0.8*M.T)/(M.T)*ones(1,Jord));
 else
  if length(M.J)<Jord-1,
   error('Order of J must be >= max(max order A, max order B)');
  end
 end
end

% Pass input and output through a non-linearity if required by model structure
if nu>0,
 x = u2x(u,M);
end
Mout    = M;
Mout.in = Mout.out;
y       = u2x(y,Mout);

%  Now Generate regressors
PHI = zeros(length(y),n+m);

%  First get regressors on X(u).
index = 1;  % Where I am up to in building up columns of PHI
for r=1:nu,  % Go through one input at a time
 for k = 0:M.nB(r), % For this input, go through all lags in numerator
  if (M.op=='d'),   
   PHI(:,index+k) = delfilter(uvec(k+1,M.nB(r)+1)',M.J,x(:,r),M.T);
  else            % Filter [p^(m_r-n)...p^(-n)]x_k
   PHI(:,index+k) = [zeros(k,1);x(1:length(x)-k,r)];   
  end
 end;
 index = index+k+1;
end;

%  Then get regressors on y.
for k = 1:n
 if (M.op=='d')
  PHI(:,m+k) = -delfilter(uvec(k+1,n+1)',M.J,y,M.T);
 else
  PHI(:,m+k) = -[zeros(k,1);y(1:length(y)-k)];
 end;
end;

%  Throw away bits where zeros were added plus extra as specified by user.
phi = PHI(OPT.n+1:end,:);
z   = y(OPT.n+1:end);

%  Can now calculate least squares estimate
if ~OPT.filt  % Don't do if this routine is only being used to generate phi
 
 switch lower(OPT.alg.type)
  case {'block'}, G=M; %Use robust inversion
   [U,S,V]=svd(phi,0); S=diag(S); r=sum(S>eps*S(1));
   G.th=V(:,1:r)*((U(:,1:r)'*z)./S(1:r));
  case {'rls','lms','kf','ct'},
   g=recur(y(:),PHI,OPT); % Solve recursively
   G=M; G.th=g.th; G.th_hist=g.th_hist; G.pe=g.pe;
  otherwise
   warning('Algorithm type not recognised, resetting to "block"');
   OPT.alg.type = 'block';  G=M; G.th = phi\z;   % Block solution if alg not recognised.
 end;
 
 %  Extract numerator and denominator polys (A & B) from parameter vector G.th
 %  Make sure we pretend there are NO non-linearities
 for i=1:nu,
  G.in(i).type = 'linear';
 end
 G.out.type = 'linear';
 G.type     = 'arx';
 if nu<1, % Special case of AR modelling
  G.B=[];
  G.nB=0;
 end
 
 %Make sure that theta has been mapped to the model
 G = theta2m(G.th,G,1);
 
 % Get Estimate of white noise variance by sample variance of residuals
 pe    = z - phi*G.th; 
 G.var = pe'*pe/length(pe);
 
 % Fill in rest of information about estimated model
 G.T     = M.T; 
 G.w     = M.w; 
 G.op    = M.op;  
 G.type  = 'arx';
 G.in    = M.in; 
 G.out   = M.out; 
 G.phi   = PHI; 
 G.delay = M.delay;
 G.T     = M.T; 
 G.w     = M.w; 
 G.op    = M.op; 
 G.delay = M.delay;
 G.type  = 'arx'; 
 G.disp.legend = 'Estimated ARX model';
 
 if ~OPT.fast  % Only do this if we have plenty of time
  
  % Get estimate of parameter covariance matrix
  G.P = G.var*pinv(phi'*phi);
  
  % Work out standard deviations on estimated parameters
  P = real(sqrt(abs(diag(G.P)))); P=P(:)';
  G.SD.A(r,:) = P(length(G.th)-n+1:length(G.th));
  if (nu>0)   % Are we only looking at an AR model?
   index = 1;
   for r=1:nu  % One T/F per i/o model
    G.SD.B(r,:)=[P(index:index+M.nB(r)),zeros(1,max(M.nB)-M.nB(r))];
    index = index+M.nB(r)+1;
    G.SD.A(r,:) = P(length(G.th)-n+1:length(G.th));
   end;
  else
   G.SD.A(1,:) = P(length(G.th)-n+1:length(G.th));
  end;
  
  % Return estimated model in ss form as well
  G = tftoss(G);
  
 end; % End of check on OPT.fast
 
 % G.alg='block'; % Record that block solution was used
 G.alg = OPT.alg;

else  % What happens if OPT.filt is set
 G.phi = PHI; G.in = M.in; G.phi = PHI;
end;

%Make sure the model type is AR when #inputs=0
if nu==0,
 G.type        = 'ar';
 G.disp.legend = 'Estimated AR model';
end

%Record that VN.m should be used for validation purposes
G.costfcn='VN';
G.OPT = OPT;