%  Running this demo's ML estimation of MIMO system via GN-based
%  algorithm and ML estimation of same MIMO system via EM algorithm.

clear; close all;  
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T    = 1;           % Sampling period in seconds
N    = 1000;        % Number of data samples
Rvar = 1e-1*eye(2); % Measurement noise

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a true linear system
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den11     = [1 1.1 0.1];
den12     = [1 2.5 1];
den21     = [1 1 0.21];
den22     = [1 1.2 0.32];
sysc      = tf({1,3; 1 1}, {den11, den12; den21, den22});
sysd      = c2d(sysc,T,'zoh');
[A,B,C,D] = ssdata(sysd); 
nx        = size(A,1);
ny        = size(C,1);
nu        = size(B,2);
delay     = [1;3];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = randn(size(B,2),N);  % The exogenous input
% Include delays specified in model structure on inputs
for r=1:2, udel(r,:) = [zeros(1,delay(r)) u(r,1:N-delay(r))]; end;
e = (sqrtm(Rvar)*randn(size(C,1),N)); % The measurement noise sequence
K = dlqe(A,eye(nx),C,eye(nx),0.1*eye(ny));
x = zeros(nx,N+1); 
y = zeros(ny,N);
for t=1:N,
	y(:,t)   = C*x(:,t) + D*udel(:,t) +   e(:,t);  %Simulate output
	x(:,t+1) = A*x(:,t) + B*udel(:,t) + K*e(:,t);  %Simulate state with innovations structure
end
Z.y = y; 
Z.u = u;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify model structure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.nx    = nx; 
M.delay = delay; 
M.estX1 = 0;
M.T     = T;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedure runs
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.dsp   = dsp; 
OPT.alg   = 'em';
OPT.miter = 100; 
opt       = OPT; 
opt.alg   = 'gn'; 
opt.cost  = 'det';   % Maximum likelihood criterion for GN search.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp, disp('Finding subspace based estimate....'); end
Gsid=sid(Z,M,OPT);
if dsp, disp('Finding ML estimate via EM....'); end
G=est(Z,M,OPT);
if dsp, disp('Finding ML estimate via gn search....'); end
Ggn=est(Z,M,opt);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,
 Gt.ss.A  = A;
 Gt.ss.B  = B;
 Gt.ss.C  = C;
 Gt.ss.D  = D;
 Gt.ss.K  = K;
 Gt.T     = T;
 Gt.w     = G.w;
 Gt.delay = delay;
 Gt.disp.colour='b';
 Gt.disp.legend = 'True Response';
 
 shownyq(Gt,G,Gsid,Ggn);	
end

echo off;

if dm
 disp('  ')
 disp('---------------------------------------------------------------------')
 disp('  ')
 disp('You now have access to the MATLAB workspace so that you may examine')
 disp('the results of this simulation.  To return to the demos, type "dbcont"')
 disp(' ')
 keyboard;
end;

