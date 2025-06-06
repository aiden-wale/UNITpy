%  Running this demo's ML estimation of grey-box parametrized MIMO model
%  structure system via GN-based algorithm.

clear; close all;  
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N    = 500;       % Number of data samples
T    = 1;         % Sampling Period
Qvar = 0;         % Variance of white state noise
Rvar = 1e-2;      % Measurement noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a linear system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den11 = [1 1.1 0.1];
den12 = [1 2.5 1];
den21 = [1 1 0.21];
den22 = [1 1.2 0.32];
sysc = tf({1,3; 1 1}, {den11, den12; den21, den22});
sysd = c2d(sysc,T,'zoh');
[A,B,C,D] = ssdata(sysd); [nx,dummy] = size(A);
sysd = ss(A,[B eye(size(A)) zeros(size(B,1),size(C,1))],C,[D zeros(size(D,1),size(A,2)) eye(size(C,1))],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t    = 0:1:N-1; 
u1   = sign(sin(2*pi*t/N)); 
u2   = sign(sin(5*pi*t/N));
u = [u1(:),u2(:)];                     % The exogenous input
w = (sqrt(Qvar)*randn(size(A,1),N))';  % State noise sequence
v = (sqrt(Rvar)*randn(size(C,1),N))';  % The measurement noise sequence
x0 = 0*ones(nx,1);
[y,t,x] = lsim(sysd,[u w v],[0:T:N*T-1],x0);
Z.y = y.'; Z.u = u.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify model structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.w     = logspace(-3,pi,300); 
M.A     = nx; 
M.type  = 'ss';
M.par   = 'struct';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedures run
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oss.dsp  = dsp;      % Subspace estimation options
oss.alg  = 'sid';

ogn.dsp  = dsp;      % ML via gn estimation options
ogn.alg  = 'gn'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Obtain an initial state%space structure estimate via subspace id.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp disp('Finding subspace based estimate....'); end;
Gsid=est(Z,M,oss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Convert this initial estimate into a canonical form
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to canonical form and now estimate via GN search
syscan=canon(Gsid.sysG,'modal');

% Set initial guess at parameters here
M.ss.A=syscan.A;
M.ss.B=syscan.B;
M.ss.C=syscan.C;
M.ss.D=syscan.D;
M.ss.K=[];
M.ss.F=[];
M.ss.G=[];
M.ss.X1=[]; %zeros(size(M.ss.A,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify which elements in structure are to be estimated % rest fixed
%  at initial values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set structure here (1's inidcate variables while 0's indicate fixed parameters)
M.ss.Ai=abs(M.ss.A)>0;
M.ss.Bi=abs(M.ss.B)>0;
M.ss.Ci=abs(M.ss.C)>0;
M.ss.Di=abs(M.ss.D)>0;
M.ss.Ki=[];
M.ss.Fi=[];
M.ss.Gi=[];
M.ss.X1i=[]; %ones(size(M.ss.X1))>0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp disp('Finding ML estimate via GN....'); end;
G=est(Z,M,ogn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp
 Mt.w=M.w; Mt.T=T; Mt.op='q';
 Mt.ss.A=A; 
 Mt.ss.B=B; 
 Mt.ss.C=C; 
 Mt.ss.D=D; 
 Mt.ss.K=[];
 Mt.type='ss'; 
 Gt=m2f(Mt);
 
 Gt.disp.colour='b';
 Gt.disp.legend = 'True Response';
 Gt.disp.aux='magonly';
 
 shownyq(Gt,Gsid,G);
end;

if dm
	disp('  ')
	disp('---------------------------------------------------------------------')
	disp('  ')
	disp('You now have access to the MATLAB workspace so that you may examine')
	disp('the results of this simulation.  To return to the demos, type "dbcont"')
	disp(' ')
	keyboard;
end;