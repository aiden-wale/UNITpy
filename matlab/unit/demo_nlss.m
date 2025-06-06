%  Running this demos estimation of a general nonlinear 
%  state space model via the maximum likelihood method
%  computed via the EM algorithm using particle smoother 
%  to compute E-step.
%
% This script demos the example profiled in the paper
% "System Identification of Nonlinear State-Space Models"
% T. Schon, A. Wills, B. Ninness
% Automatica Vol. 37, No. 1, pp. 39-49, 2011

clear; close all;  
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 100;    % Number of data samples

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a true system
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = [0.5; 25; 8; 0.05; sqrt(0.1); sqrt(0.1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = cos(1.2*[1:1:N]);
w = theta(5)*randn(1,N);
v = theta(6)*randn(1,N);
x = zeros(1,N);
y = zeros(1,N);
for k=1:N,
	x(k+1) = theta(1)*x(k) + theta(2)*(x(k)/(1+x(k)^2)) + theta(3)*u(k) + w(k);
	y(k)   = theta(4)*x(k)^2 + v(k);
end
Z.y = y;
Z.u = u;
Z   = startZ(Z);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify model structure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.type       = 'nlss';
M.nx         = 1;
M.nlss.init  = 'nlss_test_init';
M.nlss.estep = 'nlss_test_estep';
M.nlss.mstep = 'nlss_test_mstep';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedure runs
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.emit  = 50;     % The allowed number of EM iterations
OPT.miter = 100;    % The allowed number of iterations for each gradient search
OPT.pnum  = 1000;   % The number of particles used in the E-Step
OPT.lmax  = 20;     % The maximum number of bisections for each gradient search
OPT.tol   = 1e-2;   % Gradient tolerance for stopping search
OPT.dir   = 'bfgs'; % The type of search we are performing (NOTE: we don't have Hessians or approx. to them)
OPT.dsp   = dsp;    % Display (or not) the gradient search iterations
OPT.ngt   = 0;      % Perform a numerical gradient test to check the VN function for errors

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = est(Z,M,OPT);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wun = ones(1,length(G.thetait(1,:))); 
k=0:1:length(wun)-1;

%Now plot the theta trajectories

subplot(411)
h= plot(k,G.thetait(1,:),'b',k,theta(1)*wun,'r');
set(h,'linewidth',2)
ylabel('\theta_1');
xlabel('EM Iteration number');
title('Evolution of Parameter Estimates vs EM iteration number')

subplot(412)
h= plot(k,G.thetait(2,:),'b',k,theta(2)*wun,'r');
set(h,'linewidth',2)
ylabel('\theta_2')
xlabel('EM Iteration number')

subplot(413)
h= plot(k,G.thetait(3,:),'b',k,theta(3)*wun,'r');
set(h,'linewidth',2)
ylabel('\theta_3')
xlabel('EM Iteration number')

subplot(414)
h= plot(k,G.thetait(4,:),'b',k,theta(4)*wun,'r');
set(h,'linewidth',2)
ylabel('\theta_4')
xlabel('EM Iteration number')

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
