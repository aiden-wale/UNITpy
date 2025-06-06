%   EM_SUB: Subroutine called by EM.m that computes
%   expectation-maximisation steps for state space model 
%   structures.  This routine is not meant to be directly called by the
%   user - call EM.m instead which sets up all the information that needs
%   to be passed to EM_SUB.m
%
%    written by Brett Ninness,   School of EE & CS
%               Adrian Wills,    University of Newcastle
%                                Australia.

% Copyright (C) Brett Ninness, Adrian Wills

function G = em_sub(Z,M,OPT)

% Detect if gui is running
gui = 0; guih = [];
if isfield(OPT,'gui'),
 if ~isempty(OPT.gui)
  gui  = 1;         %GUI is running
  guih = OPT.gui;   %GUI handle
 end
end

% Start the display of information if desired
if OPT.dsp, 
 udisp('Cost: Negative Log-Likelihood',gui,guih);
end
sglines='----------------------------------------------------------------------------';

% Extract sizes of input and output from data matrix
[y,u,ny,nu,N] = Z2data(Z); y = y';

% Include delays specified in model structure on inputs
for r=1:nu,
 u(:,r) = [zeros(M.delay(r),1);u(1:N-M.delay(r),r)];
end
u = u';

% Need to determine what type of state space system we are dealing with and
% if all the necessary matrices are present
if findstr(M.type,'bilin'),
 bilin = 1;
else
 bilin = 0;
end

% If we don't have an input then we can't estimate all parts of a bilinear
% model and we can't estimate B and D
if numel(u)==0 || nu==0,
 bilin = 0; isBD = 0; u = zeros(0,N);
else
 isBD = 1;
end

% In what follows, a matrix H which is the expected value of
%
%   [      x(t)     ]   [      x(t)     ]'
%   [      u(t)     ]   [      u(t)     ] 
%   [kron(u(t),x(t))] * [kron(u(t),x(t))]
%   [     x(t+1)    ]   [     x(t+1)    ] 
%   [      y(t)     ]   [      y(t)     ]
%
% conditioned on the observed data will be computed in the E-step, and
% the M-step will be implemented by factorising this same matrix.
%
% Implementing the E-step will be achieved simply by callingthe Kalman 
% smoother routine with the optional flag OPT.ksem set:
OPT.ksem = 1;

% Get state dimension
nx = size(M.ss.A,1);

% Now setup indices into H based on Bilinear or not
R1i = 1 : nx+nu+bilin*nu*nx;
R2i = nx+nu+bilin*nu*nx+1 : nx+nu+bilin*nu*nx + nx + ny; 

% If initial A is unstable, then stabilize.
if max(abs(eig(M.ss.A)))>1.0,
 try,
  [XX,LL,KK] = dare(M.ss.A',M.ss.C',M.ss.Q,M.ss.R,M.ss.S);
  M.ss.A     = M.ss.A-KK'*M.ss.C;
 catch
  [XX,LL,KK] = dare(M.ss.A',M.ss.C',eye(nx),eye(size(M.ss.R)));
  M.ss.A     = M.ss.A-KK'*M.ss.C;
 end
end

% Save M into g, which will be updated at each iteration
G = M;

% Start main loop for EM method
cnt = 0;  LLold = 1e300;
while 1, 
 % E- step 
 % Implemented Kalman Smoother routine with OPT.ksem set to tell it to 
 % return H matrix described above.
 G      = ks(Z,G,OPT);
 
 % M-step - implemented by Cholesky factorisation of H 
 R      = triu(rchol(G.H));
 Gamma  = rfbs(R(R1i,R2i),R(R1i,R1i),0)';
 Pi     = (R(R2i,R2i)'*R(R2i,R2i))/N;
 
 % Extract A,B,C,D,F,G system matrices
 G.ss.A = Gamma(1:nx,1:nx);
 G.ss.C = Gamma(nx+1:nx+ny,1:nx);
 if isBD,
  G.ss.B = Gamma(1:nx,nx+1:nx+nu);
  G.ss.D = Gamma(nx+1:nx+ny,nx+1:nx+nu);
 end
 if bilin,
  G.ss.F = Gamma(1:nx,nx+nu+1:nx+nu+nu*nx);
  G.ss.G = Gamma(nx+1:nx+ny,nx+nu+1:nx+nu+nu*nx);
 end
 
 % Extract Q,S,R covariance matrices
 G.ss.Q  = Pi(1:nx,1:nx);
 G.ss.S  = Pi(1:nx,nx+1:nx+ny);
 G.ss.R  = Pi(nx+1:nx+ny,nx+1:nx+ny);
 
 % Extract initial state estimate X1 with cov P1
 G.ss.X1 = G.ss.xs(:,1);
 G.ss.P1 = G.ss.Ps(:,:,1)'*G.ss.Ps(:,:,1);
 
 % Update the iteration counter
 cnt = cnt + 1;
 
 % Print iteration information to the screen
 if OPT.dsp,
  if cnt==1,
   udisp(sglines,gui,guih);
   str1=sprintf('%s%13s','Iter#','Cost');
   udisp(str1,gui,guih);
  end
  udisp(sprintf('%5i%13.3e',cnt,G.LL),gui,guih);
 end
 
 % Check the stopping criteria
 if abs(G.LL-LLold) < OPT.mdec*(abs(G.LL)+abs(LLold)), 
  G.whystop = 'Termination due to Log-likelihood difference less than OPT.mdec';
  if OPT.dsp,
   udisp(sglines,gui,guih);
   udisp(G.whystop,gui,guih);
   udisp(sglines,gui,guih);
  end
  break; 
 end
 if cnt >= OPT.miter,
  G.whystop = 'Termination due to number of iterations exceeding OPT.miter';
  if OPT.dsp,
   udisp(sglines,gui,guih);
   udisp(G.whystop,gui,guih);
   udisp(sglines,gui,guih);
  end
  break;
 end
 LLold = G.LL; % Store the old log-likelihood
 
 
 % Store parts of G that we want and remove the rest
 gg=M;
 gg.ss.A=G.ss.A;
 gg.ss.C=G.ss.C;
 gg.ss.Q=G.ss.Q;
 gg.ss.S=G.ss.S;
 gg.ss.R=G.ss.R;
 gg.ss.P1=G.ss.P1;
 gg.ss.X1=G.ss.X1;
 if isBD,
  gg.ss.B = G.ss.B;
  gg.ss.D = G.ss.D;
 end
 if bilin,
  gg.ss.F = G.ss.F;
  gg.ss.G = G.ss.G;
 end
 
 clear G;
 G=gg;
end
