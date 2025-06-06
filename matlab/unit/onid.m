%   This function does least squares fitting with respect to orthogonal
%   basis vectors for observed input-output data u and y.  The orthogonal
%   basis vectors are formed from filtered versions of u.  A special case
%   of FIR modelling occurs if all the fixed poles are specified as zero.
%
%   Usage is:   G = onid(Z,M,OPT)
%
%   Where:
%
%   Z         Input output data in the form Z = [y,u] where y is a column
%             vector of output measurements and u is a matrix whose
%             columns are the input measurements - this means that MISO
%             models are catered for, for MIMO the user should conduct
%             multiple MISO estimation runs (one for each output).
%   M         Data structure which defines the model structure which
%             is to be estimated from the data as follows:
%   M.poles   Poles to be used in basis vectors.
%  	          If you include a complex pole, then don't include
%  	          its complex conjugate, since the conjugate is
%  	          automatically included for you.
%   M.delay   Number of of samples of delay to include.  In the
%             case of a MISO system, this should be a vector of delays,
%             one for each input being considered.
%   M.T       sampling period in s. (Ignored for q case) Default = 1;
%   M.w       vector of frequencies at which to calculate frequency
%             response of estimated model.  Specify in real frequency,
%             not normalised.  Default is 3 decades up to folding freq.
%   OPT.n     Number of initial samples to discard to get rid of
%  	      transients in regressor filters.
%   OPT.fast  When set to 1, then makes algorithm run maximally fast by
%             avoiding the calculation of error bounds and avoiding
%             calculation of G.GAMMA.  Default is 0.
%   OPT.filt  When set to 1, then algorithm only does filtering to
%             generate G.PHI - no least squares estimation is done.
%   OPT.alg   This structure sets the algorithm used to find the estimate
%             and also sets any parameters that are associated with
%             the algorithm. The components of this structure are as follows:
%   OPT.alg.type - determines the algorithm type used.  It may be set as:
%            'block' - the solution for the least squares
%                      estimate is found in block form (default).
%            'rls'   - the solution is found recursively
%                      the via recursive least squares algorithm.
%            'ct'    - the solution is found recursively via recursive
%                      least squares with a contant trace forced
%                      on the covariance matrix.
%            'lms'   - the solution is found recursively via
%                      least-mean-square (gradient descent) algorithm.
%            'kf'    - the solution is found recursively via
%                      a Kalman Filtering algorithm.
%   OPT.alg.P      - Initialisation of covariance matrix, default is P=10*I
%   OPT.alg.mu     - LMS gain, default is mu=0.001;
%   OPT.alg.lambda - RLS alg. `forgetting factor'. Default is lambda = 1;
%   OPT.alg.R      - Measurement noise variance used by kf alg. Default = 1;
%   OPT.alg.Q      - Parameter variance used by kf alg. Default = 0.1*I;
%   OPT.alg.th     - Initial parameter vector value.  Default = 0.
%
%   Output variables available at the end of this macro are:
%
%   G.B/G.A   estimated transfer function for model of dynamics.
%   G.G       Frequency response of estimated model for dynamics.
%   G.Ge      Matrix specifying 95% confidence regions for estimated
%             frequency response Ghat.  They may be plotted by using either
%             of the commands `shownyq(G)' or `showbode(G)'.
%   G.P       Covariance Matrix for Estimated Parameters.
%   G.phi     Matrix of regressors generated from data such that y = PHI*y
%   G.th      Model estimates as a column vector
%   G.TH      Same as above, but for MISO systems G.TH is G.th
%             re-arranged as a matrix, each column being the parameters
%             for a corresponding input.
%   G.GAMMA   Vector such that estimated model frequency
%             response = G.GAMMA*G.th;
%
%   If a recursive solution algorithm is used then also available is:
%
%   G.th_hist History of how parameter space estimates evolved.
%   G.pe      Time history of how prediction error evolved.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness.

function G = onid(Z,M,OPT)

%  Extract output and inputs along with their dimensions
Z=startZ(Z);  [y,u,ny,nu,Ny] = Z2data(Z);

% Unspecified parts of OPT -> defaults
OPT = startOPT(OPT);
if ~isfield(OPT.alg,'type'),
 OPT=rmfield(OPT,'alg');
 OPT.alg.type = 'block';
end

% Unspecified parts of M -> defaults
if (OPT.n>=Ny) 
 error('Cannot have OPT.n larger than height of Z!'); 
end
if ~exist('M','var'), 
 error('Need to specify initial model structure M!');
else
 M = startM(Z,M);
end

% Form rich test input - used for conversion to SS/TF form.
Nimp = (6*length(M.poles)+2)*(nu+1)*2; 
imp  = randn(Nimp,nu);

% Include delays specified in model structure on inputs
for r=1:nu
 u(:,r) = [zeros(M.delay(r),1);u(1:Ny-M.delay(r),r)];
 % imp(:,r) = [zeros(M.delay(r),1);imp(1:Nimp-M.delay(r),r)];
end;

% Pass input and output through a non-linearity if required by model structure
x = u2x(u,M); %Mout=M; Mout.in=Mout.out; y = u2x(y,Mout); %WTF? Why pass output through NL

%Need to make sure that the user has not supplied the complex conjugate
idx = [];
for ii=1:length(M.poles),
 if ~isreal(M.poles(ii)),
  for jj=ii+1:length(M.poles),
   if abs(M.poles(ii)-conj(M.poles(jj)))<eps,
    idx = [idx jj];
   end
  end
 end
end
M.poles(idx) = [];

%  Now check to see if an FIR model has been specified in which case we can
%  save ourselves a lot of work
if ( max(abs(M.poles)) > 0 )
 
 %  Now generate the regressors as specified by the poles.  Also generate
 %  the matrix GAMMA which allows the frequency response of the model
 %  to be easily calculated.
 PHI = zeros(Ny,(length(M.poles)+length(find(imag(M.poles))))*nu);
 if ~OPT.fast GAMMA = zeros(length(M.w),length(M.poles) + length(find(imag(M.poles)))); end;
 ww = exp(j*M.w*M.T);
 xap = x; ip=imp;       %  All-pass filtered version of input and impules
 gap = ones(size(ww));  %  Freq response of all-pass filter
 pindex = 1;            %  What column of phi we are up to.
 
 for k = 1:length(M.poles)
  %  First test to see if the pole is complex or not.
  if ( abs( imag( M.poles(k) ) ) > 0 )
   %  If a pole is complex we have to include its conjugate
   d1 = real(conv([1,-M.poles(k)],[1,-M.poles(k)']));
   alpha = 2*real(M.poles(k))/(1+abs(M.poles(k))^2);
   %  Go for conventional `2-parameter' Kautz Numerator as default.
   num1=sqrt((1-alpha^2)*(1+abs(M.poles(k))^2)*(1-abs(M.poles(k))^2))*[0,0,1];
   num2=sqrt((1+abs(M.poles(k))^2)*(1-abs(M.poles(k))^2))*[0,1,-alpha];
   phi = filter(num1,d1,xap);
   PHI(:,(pindex-1)*nu+1:(pindex-1)*nu+nu) = phi;
   PSI(:,(pindex-1)*nu+1:(pindex-1)*nu+nu) = filter(num1,d1,ip);
   if ~OPT.fast GAMMA(:,pindex) = polyval(num1,ww(:))./polyval(d1,ww(:)).*gap(:); end;
   pindex = pindex+1;
   phi = filter(num2,d1,xap);
   PHI(:,(pindex-1)*nu+1:(pindex-1)*nu+nu) = phi;
   PSI(:,(pindex-1)*nu+1:(pindex-1)*nu+nu) = filter(num2,d1,ip);
   if ~OPT.fast
    GAMMA(:,pindex) = polyval(num2,ww(:))./polyval(d1,ww(:)).*gap(:);
    gap = (polyval([d1(3:-1:1)'],ww(:))./polyval(d1,ww(:))).*gap(:);
   end;
   pindex = pindex+1;
   xap=filter([d1(3:-1:1)'],d1,xap);  ip=filter([d1(3:-1:1)'],d1,ip);
  else
   n1 = [0,sqrt(1-abs(M.poles(k))^2)];
   d1 = [1,-M.poles(k)];
   phi = filter(n1,d1,xap);
   PHI(:,(pindex-1)*nu+1:(pindex-1)*nu+nu) = phi;
   PSI(:,(pindex-1)*nu+1:(pindex-1)*nu+nu) = filter(n1,d1,ip);
   if ~OPT.fast
    GAMMA(:,pindex) = (polyval(n1,ww(:))./polyval(d1,ww(:))).*gap(:);
    gap = (polyval(fliplr(d1),ww(:))./polyval(d1,ww(:))).*gap(:);
   end;
   pindex = pindex+1;
   xap=filter(fliplr(d1),d1,xap);
   ip=filter(fliplr(d1),d1,ip);
  end;  %  End of check for complex pole or not
  
 end;  %  End of iteration through poles
 
else  %  Special case of FIR
 r = length(M.poles)+1; PHI = zeros(Ny,r*nu);
 [rr,tt] = size(imp); PSI = zeros(rr,r*nu);
 for t=1:nu
  index = [zeros(1,t-1),1,zeros(1,nu-t)];
  phi = toeplitz(x(:,t),[x(1,t),zeros(1,length(M.poles))]);
  psi = toeplitz(imp(:,t),[imp(1,t),zeros(1,length(M.poles))]);
  PHI = PHI + kron(phi,index); 
  PSI = PSI + kron(psi,index);
 end;
 if ~OPT.fast GAMMA=exp(-j*kron(M.w(:)*M.T,0:1:r-1)); end;
end;

%  Now get the estimate via least squares (block) or some recursive method
if ~OPT.filt
 switch lower(OPT.alg.type)
  case {'block'}, G.th = PHI(OPT.n+1:length(y),:)\y(OPT.n+1:length(y)); % Block solution
  case {'rls','lms','kf','ct'},    G = recur(y(:),PHI,OPT);             % Solve recursively
  otherwise
   warning('Algorithm type not recognised, resetting to "block"')
   OPT.alg.type = 'block';
   G.th = PHI(OPT.n+1:length(y),:)\y(OPT.n+1:length(y)); % Block solution if alg not recognised.
 end;   % End of switch statement
end;    % End of check on OPT.filt

% Load up output with model properties
G.phi   = PHI; 
G.poles = M.poles; 
G.type  = 'fir'; 
G.op    = 'q';
G.T     = M.T; 
G.w     = M.w; 
G.in    = M.in; 
G.out   = M.out; 
G.nA    = M.nA; 
G.nB    = M.nA-1; 
G.delay = M.delay;
G.nu    = M.nu; 
G.ny    = M.ny;

% Add legend for prospective plotting
G.disp.legend=['Estimated (Generalised) FIR model'];

G.alg='block'; % Record that block solution was used

% Do luxurious extras (if not doing fast version)

if ~OPT.fast
 
 % Get estimated model in state space form:
 m.nx   = length(M.poles)+length(find(imag(M.poles)));
 m.op   = 'q'; 
 m.T    = 1; 
 OO.horizon = m.nx+2;
 yimp   = PSI*G.th; 
 g      = sid([yimp(:),imp],m,OO);
 G.ss.A = g.ss.A; 
 G.ss.B = g.ss.B; 
 G.ss.C = g.ss.C; 
 G.ss.D = g.ss.D;
 pdel   = exp((-j*M.w(:)*M.T)*G.delay');
 
 for r=1:nu  % One freq response and T/F per i/o model
  index      = r:nu:(length(M.poles)+length(find(imag(M.poles))))*nu;
  G.G(1,r,:) = GAMMA*G.th(index); 
  pp(1,1,:)  = pdel(:,r); 
  G.G(1,r,:) = G.G(1,r,:).*pp;
  [G.B(r,:),G.A(r,:)] = ss2tf(g.ss.A,g.ss.B(:,r),g.ss.C,g.ss.D(1,r),1);
 end;
 
 % Parameter space variance of estimates:
 G.C = []; G.D = []; G.GAMMA = GAMMA;
 pe = y(OPT.n+1:length(y)) - PHI(OPT.n+1:length(y),:)*G.th;
 G.var = pe'*pe/length(pe);
 G.P = G.var*pinv(PHI(OPT.n+1:length(y),:)'*PHI(OPT.n+1:length(y),:));
 
 % Now load up matrix specifying standard deviations
 P = real(sqrt(diag(G.P))); P = P(:); d = length(M.poles)+length(find(imag(M.poles)));
 for r=1:nu G.SD.th(:,r) = P((r-1)*d+1:r*d); G.TH(:,r) = G.th((r-1)*d+1:r*d); end;
end;


%Record that validate should use VN as the cost function to obtain
%prediction errors
G.costfcn = 'VN';

G.OPT = OPT;



