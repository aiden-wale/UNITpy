%  Running this tests estimation of randomly chosen 8'th order 
%  Bilinear system via both PEM criterion with gradient based 
%  search, and ML criterion using the EM algorithm

clear; close all; global dm; if isempty(dm), clear global dm; dm=0; end
global trans;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N    = 500;    % Number of data samples
Qvar = 1e-2;   % State noise variance
Rvar = 1e-1;   % Measurement noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Generate randomly chosen bilinear system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=8;     % State dimension
m=2;     % Number of inputs
p=2;     % Number of outputs
g=drss(n,p,m);
[A,B,C,D]=ssdata(g); A=0.99*A;
F=0.1*randn(n,n*m);
G=0.1*randn(p,n*m);
Q=sqrt(Qvar)*eye(n); R = sqrt(Rvar)*eye(p); S = zeros(n,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=randn(m,N); 
w=Q*randn(n,N); v=R*randn(p,N);
x=zeros(n,N); y=zeros(p,N);
for k=1:N,
 ukx = kron(u(:,k),x(:,k));
 x(:,k+1) = A*x(:,k) + B*u(:,k) + F*ukx + w(:,k);
 y(:,k) = C*x(:,k) + D*u(:,k) + G*ukx + v(:,k);
end
Z.y = y; Z.u = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure as random perturbation from truth                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.A=n; M.type='bilinear';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how 
%  estimation procedure runs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Oem.dsp=1; Oem.miter=500; Oem.alg='em';
Ogn.dsp=1; Ogn.miter=500; Ogn.alg='gn';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gem=est(Z,M,Oem);
Ggn=est(Z,M,Ogn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display results in - validation on observed data set
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

validate(Z,Gem);
%disp('Current plots show EM results.'); disp('Press enter for GN results');pause;
validate(Z,Ggn);

if dm
 disp('  ')
 disp('---------------------------------------------------------------------')
 disp('  ')
 disp('You now have access to the MATLAB workspace so that you may examine')
 disp('the results of this simulation.  To return to the demos, type "dbcont"')
 disp(' ')
 keyboard; 
end;
