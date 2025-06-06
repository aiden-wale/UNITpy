%  Running this demos EM estimation of a randomly chosen multivariable
%  system using subspace-based estimation, least squares estimation via
%  Gauss-Newton search, and ML estimation via EM-algorithm search.

clear; close all; global dm; if isempty(dm), clear global dm; dm=0; end
global trans;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions                    
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N    = 200;    % Number of data samples
T    = 1;      % Sampling Period
Qvar = 1e-2;   % Variance of white state noise 
Rvar = 1e-0;   % Variance of white measurement noise 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a linear system
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a random system
ny = 2; % Number of system outputs 
nu = 2; % Number of system inputs
nx = 5; % Number of system states/ Model order

% Generate a stable discrete-time state space system in innovations form
sys=drss(nx,ny,nu); [A,B,C,D]=ssdata(sys); A=0.99*A;
Q=1*eye(nx); R=1*eye(ny); [P,L,K]=dare(A',C',Q,R);  K=K';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sdv = diag(1e-1*ones(ny,1));
u   = randn(nu,N);
v   = sdv*randn(ny,N);
R   = sdv*sdv;

%Use lsim to generate output

sys = ss([A-K*C],[B-K*D -K],C,[D zeros(ny)],-1); 
y   = lsim(sys,[u;v])'+v;
Z.y = y; Z.u=u;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedure runs
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oss.dsp = 1;
oss.alg = 'sid';

oem.dsp = 1; 
oem.alg = 'em'; 

ogn.dsp  = 1;
ogn.alg  = 'gn'; 
opt.cost = 'trace';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init = 2; % 1=random, 2=subspace
switch init, 
 case 1,
  % Perturb initial system
  alp=0.01; beta=1; gamma=0;
  M.ss.A=beta*(gamma*A + alp*rand(size(A))); 
  M.ss.B=beta*(gamma*B + alp*rand(size(B))); 
  M.ss.C=beta*(gamma*C + alp*rand(size(C)));  
  M.ss.D=beta*(gamma*D + alp*rand(size(D)));  
  M.ss.K=beta*(gamma*K + alp*rand(size(K)));  
  M.ss.Q=100*eye(nx); M.ss.R=0.01*eye(ny);
 case 2,
  %Use LTI subspace method
  clear M;
  M.A = nx;
  O.alg='n4sid';
  M = sid(Z,M,O);
  M.ss.Q=1*eye(nx); M.ss.R=0.01*eye(ny); M.ss.S=0*M.ss.S;
end

M.w = logspace(-3,pi,10000); 
M.type='ss';

disp('Finding subspace based estimate....')
gss=est(Z,M,oss);
disp('Done');

disp('Finding ML estimate via EM....')
gem=est(Z,gss,oem);
disp('Done');
disp('Finding ML estimate via gn search....')
ggn=est(Z,gss,ogn);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mt.w=M.w; Mt.delay=ggn.delay; Mt.T=T; Mt.op='q';
Mt.ss.A=A-K*C; Mt.ss.B=B-K*D; Mt.ss.C=C; Mt.ss.D=D; Mt.ss.K=[];
Mt.type='ss'; 
Gt=m2f(Mt);

Gt.w = M.w; Gt.T=ggn.T;  Gt.disp.colour='b'; 
Gt.disp.legend = 'True Response';
Gt.disp.aux='magonly';
Gs.disp.legen='slow em';
Gs.disp.linestyle=':';

shownyq(Gt,gem,gss,ggn);

if dm
 disp('  ')
 disp('---------------------------------------------------------------------')
 disp('  ')
 disp('You now have access to the MATLAB workspace so that you may examine')
 disp('the results of this simulation.  To return to the demos, type "dbcont"')
 disp(' ')
 keyboard; 
end;




