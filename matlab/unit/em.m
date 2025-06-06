%   EM: Function to find maximum likelihood estimates of state space mdoel
%   by means of the Expectation-Maximisation (EM) algorithm.  This
%   function is not meant to be called directly.  Rather, it is intended
%   to be called by using the est.m algorithm, with OPT.alg='em' specified.
%
%   Usage is:  G = em(Z,M,OPT);
%
%   where
%
%     Z:        Input-Output data in one of two forms.  The standard form
%               is for it to be a record with elements Z.y and Z.u, each
%               of which are matrices with number of rows equal to the
%               number of data samples, and number of columns equal (respectively)
%               to the number of outputs and the number of inputs.  On
%               the other hand, Z can be a matrix of the form Z = [y,u]
%               where it is assumed that y is a column vector of output
%               measurements and u is a matrix whose columns are the
%               input measurements; in this latter MISO models are
%               being considered.
% M.ss.A,B,C:   Initial state space model structure guess.  These *cannot*
%      D,F,G    be specified simply as integer orders, they must be an
%               actual full state space system parameterisation.  If
%               M.ss.D is set as an empty matrix, then a feedthrough term
%               is not estimated, andthe returned estimate G.ss.D is set
%               at zero;
%     M.delay:  Number of samples of delay to include. In the
%               case of a MIMO system, this should be a vector of delays,
%               one for each input being considered.
%     M.op:     set to 'q' for shift and 'd' for delta.  Default = 'q'.
%     M.T:      sampling period (ignored for q operator case).  Default=1
%     M.w:      vector of frequencies at which to calculate frequency
%               response of estimated model.  Specify in real frequency,
%               not normalised.  Default is 3 decades up to folding freq.
%     M.type:   If set to 'ss' (default) a linear state space model structure
%               is estimated (M.ss.F, M.ss.G then irrelevant).  If set to
%               'bilinear', then a bilinear state space structure is estimated.
%    OPT:       Data structure which defines options for the estimation
%               algorithm as follows:
%    OPT.miter: Maximum number of iterations in search for minimum.  Default = 30.
%    OPT.dsp:   Control of output to screen 0=>quiet,1=>verbose.  Default = 0
%    OPT.mdec:  Minimum relative decrease of cost before search is
%               terminated.  Default = 1e-4.
%    G:         Data structure which specifies the estimated model as
%               follows:
%  G.A, G.B     Matrices definining the estimated transfer function model.
%  G.C, G.D     For SISO systems, these element are row vectors defining
%               co-efficients of increasing powers of M.op^-1.  For MISO,
%               they are matrices of rows, the k't row pertaining to the
%               k'th input.  For MIMO, they are 3 dim matrices with the
%               row (k,:,m) defining the transfer function from the k'th
%               input to the m'th output.
% G.ss.A,B,C:   [A,B,C,D,F,G] matrices/vectors defining estimated state space
%       F,G:    model.
% G.ss.X1, P1:  Estimate of initial conditions and covariance of that estimate.
%    G.G:       Matrix of frequency responses.  If the system has multiple
%               inputs and multpile outputs, then this matrix is 3
%               dimensional, with one `page' per output, and the i'th
%               column of the j'th page (ie G.G(i,j,:)) being the
%               frequency response from input i to ouput j.
%    G.mse:     Evolution of mean square cost decrease;
%    G.LL:      Evolution of log likelihood increase.
%
%    written by Brett Ninness,   School of EE & CS
%               Adrian Wills,    University of Newcastle
%                                Australia.

% Copyright (C) Brett Ninness, Adrian Wills

function G = em(Z,M,OPT);

% Unspecified parts of OPT -> defaults
if ~exist('OPT'), 
 OPT = startOPT([]); 
else
 OPT = startOPT(OPT); 
end

% Copy EM algorithm result into output structure G
G=M; % Input specs on model type etc are copied to output

%switch based on data type
switch Z.type,

 case 'time',
  % Extract sizes of input and output from data matrix
  [y,u,ny,nu,N] = Z2data(Z);

  % Include delays specified in model structure on inputs
  for r=1:nu,
   u(:,r) = [zeros(M.delay(r),1);u(1:N-M.delay(r),r)];
  end
  Zem.u=u';
  Zem.y=y';

  % Check that initialisations make sense and convert to form suitable for mex file
  m=M.ss; m.type=M.type; m.in=M.in; [m,nx] = sschk(m,nu,ny);

  % Check initialisation of covariance matrices makes sense
  while 1,
   Pi = [m.Q m.S;m.S' m.R];
   if any(eig(Pi)<eps), % A little ad-hoc? Check for non-pos [Q,S;S^T R] and add seasoning if necessary
    m.Q=m.Q+max(1e-5,norm(m.Q))*eye(size(m.Q));
    m.R=m.R+max(1e-5,norm(m.R))*eye(size(m.R));
   else
    break;
   end
  end
  m.Q=Pi(1:nx,1:nx); m.S=Pi(1:nx,nx+1:end); m.R=Pi(nx+1:end,nx+1:end);
  M.ss.Q=m.Q; M.ss.S=m.S; M.ss.R=m.R;
  
  if OPT.dsp,
   disp('Algorithm: Expectation-Maximisation');
  end

  % Call the EM routine implemented as a mex file for speed
  if strcmp(M.type,'ss')
   % Check to see if non-linearity on any input
   lin=1; for k=1:length(M.in) lin=lin*strcmp(lower(M.in(k).type),'linear'); end;
   if ~lin  % OK, one input (at least) was non-linear => Hammerstein version
    g=em_hamm(Zem,m,OPT);
    % State/measurement noise cov + ic estimates -> G
    G.ss.Q=g.ss.Q;G.ss.R=g.ss.R;G.ss.S=g.ss.S;G.ss.P1=g.ss.P1;G.ss.X1=g.ss.X1;
    %Record that validate should use VN as cost function
    G.costfcn = 'VN';
   else
    g=em_sub(Z,M,OPT);
    G=g;
    G.costfcn = 'VNss';
   end;
  elseif any(strcmp(M.type,{'bilin','bilinear'}));
   g=em_sub(Z,M,OPT);
   G=g;
   G.costfcn = 'VNss';
  elseif any(strcmp(M.type,{'nbj'}))
   g=em_hamm(Z,M,OPT);
   % State/measurement noise cov + ic estimates -> G
   G.ss.Q=g.ss.Q;G.ss.R=g.ss.R;G.ss.S=g.ss.S;G.ss.P1=g.ss.P1;G.ss.X1=g.ss.X1;
   %Record that validate should use VN as cost function
   G.costfcn = 'VN';
   % Estimates of system matrices -> G
   G.ss.A=g.A;G.ss.B=g.B;G.ss.C=g.C;G.ss.D=g.D;
   G.ss.F=[]; G.ss.G=[];
  else
   error('M.type is unknown');
  end

  % Evolution of Likelihood and mean square error -> G
  G.mse=g.PE; G.LL=g.LL;

   % Compute estimated innovations variance -> G
  [P,L,K] = dare(G.ss.A',G.ss.C',G.ss.Q,G.ss.R,G.ss.S);
  G.var   = G.ss.C*P*G.ss.C'+G.ss.R;
  G.ss.K  = K';
  
  
 case 'frequency',
  [y,w,ny,nu,N] = Z2data(Z);  
  G = fem(Z,M,OPT);
  
end
  
% Convert linear components from state space to transfer function form
G = sstotf(G);

% Add legend for prospective plotting
G.disp.legend=['Estimated ',upper(G.type),' model via EM'];

% Record that EM algorithm was used
G.alg='em';

% Pass out any estimate of input non-linearity
if nu>0
 G.in = M.in;
 %  for k=1:nu, G.in(k).type='linear'; G.in(k).neta=0; end;
else  % The ARMA case
 G.in.type='linear'; G.in.neta=0;
end;
if ~isfield(M,'out'), 
 G.out.type='linear'; 
 G.out.neta=0; 
end
