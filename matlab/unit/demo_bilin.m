%  Running this demos estimation of multivariable bilinear system
%  via both PEM criterion with gradient based search, and ML 
%  criterion using the EM algorithm

clear; close all;
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N    = 1000;    % Number of data samples
Rvar = 1e-4;   % Measurement noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Generate Favoreel bilinear test system
%  Ref: IEEE TAC 44(6), pp 1157-1165, 1999
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=[0.5,0;0, 0.3]; B=[0,1;-1,0]; C=[1,0;0,2]; D=eye(2,2);
F=[0.6,0,0.2,0;0,0.4,0,0.5]; H= zeros(2,4);
R = sqrt(Rvar)*eye(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u      = randn(2,N);
v      = R*randn(2,N);
x      = zeros(2,N);
y      = zeros(2,N);
x(:,1) = 0*0.5*[1;-1];
for k=1:N,
	ukx = kron(u(:,k),x(:,k));
	x(:,k+1) = A*x(:,k) + B*u(:,k) + F*ukx;
	y(:,k) = C*x(:,k) + D*u(:,k) + H*ukx + v(:,k);
end
Z.y = y; 
Z.u = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.nx   = 2;             % Order of state dimension
M.type = 'bilinear';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how
%  estimation procedure runs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Oem.dsp=dsp; 
Oem.miter=200; 
Oem.alg='em';

Ogn.dsp=dsp; 
Ogn.miter=200; 
Ogn.alg='gn';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gem=est(Z,M,Oem);
Ggn=est(Z,M,Ogn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display validation results on observed data set
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dsp,
	validate(Z,Gem);
	validate(Z,Ggn);
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
