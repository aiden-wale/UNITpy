%  Running this demos the estimation of a multivariable 
%  state space time series model using the Maximum Likelihood
%  method computed by the EM algorith and gradient based search

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

N   = 10000;     % Number of data samples
T   = 1;        % Sampling Period
var = 0.1;      % Innovations variance

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a true (randomly chosen) linear system
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ny = 2; % Number of system outputs
nx = 2; % Number of system states

% Generate a randomly chosen stable discrete-time state space system in innovations form
sys       = drss(nx,ny,ny); 
[A,K,C,D] = ssdata(sys);

% This gives stable poles chosen randomly, but might have non-minimum
% phase zeros.  Get a valid spectral factor via dare solution with big Q.

P = dare(A',C',eye(nx,nx),eye(ny,ny)); 
K = (A*P*C')*inv(C*P*C'+eye(ny,ny))';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e   = sqrt(var)*randn(N,ny);  % Innovations
x   = ltitr(A,K,e);
Z.y = C*x' + e'; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify model structure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.nx   = nx; 
M.type ='ss'; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedure runs
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.dsp   = dsp; 
OPT.alg   = 'em';
OPT.fast  = 1;
opt       = OPT; 
opt.alg   = 'gn'; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,	disp('Finding subspace based estimate....'); end
Gsid=sid(Z,M,OPT);
if dsp, disp('Done'); disp('Finding ML estimate via EM....'); end
G=est(Z,Gsid,OPT);
if dsp, disp('Done'); disp('Finding ML estimate via gn search....'); end
Ggn=est(Z,Gsid,opt);
if dsp, disp('Done'); end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dsp,
	Gt.ss.A = A;
	Gt.ss.C = C;
	Gt.ss.D = D;
	Gt.ss.K = K;
    Gt.w     = G.w;    
	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';

	showbode(Gt,Gsid,Ggn,G);
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