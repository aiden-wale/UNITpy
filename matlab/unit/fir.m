%   This function does least squares fitting for FIR models based on 
%   observed input-output data u and y.
%
%   Usage is:   G = fir(Z,M,OPT)
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
%   M.delay   Number of of samples of delay to include.  In the
%             case of a MISO system, this should be a vector of delays,
%             one for each input being considered.
%   M.T       sampling period in s. (Ignored for q case) Default = 1;
%   M.w       vector of frequencies at which to calculate frequency
%             response of estimated model.  Specify in real frequency,
%             not normalised.  Default is 3 decades up to folding freq.
%   OPT.n     Number of initial samples to discard.
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
%   G.B       estimated transfer function for model of dynamics.
%   G.P       Covariance Matrix for Estimated Parameters.
%   G.th      Model estimates as a column vector
%   G.phi     Matrix of regressors generated from data such that y =
%             PHI*G.th

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

function G = fir(Z,M,OPT)

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

% Include delays specified in model structure on inputs
for r=1:nu
 u(:,r) = [zeros(M.delay(r),1);u(1:Ny-M.delay(r),r)];
end;

% Pass input through a non-linearity if required by model structure
x = u2x(u,M);

% For regressor matrix
PHI = zeros(Ny,sum(M.nB+1));
idx = 0;
for idu = 1:nu
 PHI(:,idx+1:idx+M.nB(idu)+1) = toeplitz(x(:,idu),[x(1,idu),zeros(1,M.nB(idu))]);
 idx = idx + M.nB(idu) + 1;
end;

% Save initial model into G
G = M;

%  Now get the estimate via least squares (block) or some recursive method
if ~OPT.filt
 switch lower(OPT.alg.type)
  case {'block'}, 
   G.th = PHI(OPT.n+1:length(y),:)\y(OPT.n+1:length(y)); % Block solution
  case {'rls','lms','kf','ct'},    
   g    = recur(y(:),PHI,OPT);             % Solve recursively
   fn   = fieldnames(g);
   for i=1:length(fn),
    if ~isfield(G,fn{i}),
     G = setfield(G,fn{i},getfield(g,fn{i}));
    elseif isempty(getfield(G,fn{i}))
     G = setfield(G,fn{i},getfield(g,fn{i}));
    end
   end
   
  otherwise
   warning('Algorithm type not recognised, resetting to "block"')
   OPT.alg.type = 'block';
   G.th = PHI(OPT.n+1:length(y),:)\y(OPT.n+1:length(y)); % Block solution if alg not recognised.
 end;   % End of switch statement
end;    % End of check on OPT.filt

% Put theta into model structure
idx = 0;
mxB = max(M.nB);
G.B = zeros(nu,mxB+1);
for idu = 1:nu
 G.B(idu,:) = [G.th(idx+1:idx+M.nB(idu)+1)' zeros(1,mxB-M.nB(idu))];
 idx        = idx + M.nB(idu) + 1;
end;


% Do luxurious extras (if not doing fast version)
if ~OPT.fast
 % Parameter space variance of estimates:
 pe    = y(OPT.n+1:length(y)) - PHI(OPT.n+1:length(y),:)*G.th;
 G.var = pe'*pe/length(pe);
 G.P   = G.var*pinv(PHI(OPT.n+1:length(y),:)'*PHI(OPT.n+1:length(y),:));
 
 % Now load up matrix specifying standard deviations
 P = real(sqrt(diag(G.P))); 
 P = P(:); 
 d = 0;
 for r=1:nu 
  G.SD.th(:,r) = [P(d+1:d+M.nB(r)+1);  zeros(mxB-M.nB(r),1)];
  d = d + M.nB(r)+1;
 end
end

% Load up output with model properties
G.phi   = PHI; 

% Record that validate should use VN as the cost function to obtain
% prediction errors
G.costfcn = 'VN';
G.OPT = OPT;

% Add legend for prospective plotting
G.disp.legend=['Estimated n_b=',int2str(max(M.nB)),'th order FIR model'];
G.alg='block'; % Record that block solution was used


