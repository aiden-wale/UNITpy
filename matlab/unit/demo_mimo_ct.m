%  Running this demo's PEM estimation of continuous time state space
%  model via gradient-based search.

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
 
T  = 1e-3;                % Sampling period for simulation 
                          % of true system (sec)
Tf = 256;                 % Duration of data record (sec)
Ns = round(Tf/T);         % Implied number of samples
ts = linspace(0,Tf-T,Ns); % Sampel time points

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a true linear system
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mt.ss.A  = [-0.5 1; -0.6 0];
Mt.ss.B  = [1;1];
Mt.ss.C  = [1 0];
Mt.ss.D  = [];
K = 1e-3*[1;1];
Mt.ss.K  = K;
Mt.ss.Q  = K*K';
Mt.ss.S  = K;
Mt.ss.R  = 1;
Mt.ss.X1 = zeros(size(Mt.ss.A,1),1);
Mt.T = T;
Mt.op = 's';

% Obtain its sampled data equivalent assuming zero order hold on input
% and integrated sampling on the output

Ms = samplek(Mt);
A=Ms.ss.A; B=Ms.ss.B; C=Ms.ss.C; D=Ms.ss.D;
Q=Ms.ss.Q; S=Ms.ss.S; R=Ms.ss.R;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

regular = 1;   % True system simulated at T second time scale
               % can be either regularly or irregularly sampled.

% First, generate the input and sampling indices

if regular,    % Regular sampling requested       

 Ts  = 0.2;          % Regular sampling interval (sec)
 s   = round(Ts/T);  % Factor by which simulation sampling is faster
 N   = floor(Tf/Ts); % Number of samples available for estimation
 idx = 1:s:Ns;       % Indexes for regular sampling of true system
 
 % Generate the input (regulararly spaced samples)
 u  = randn(1,N); 
 tt = 0:1:N-1; u = sign(sin(4*pi*tt/N));
 
 uc = u(kron([1:N],ones(1,s)));  % Upsample the input to correct the time base

else   % Irregular sampling requested

 Ts  = 2;
 s   = round(Ts/T);  % Number of times greater to sample discrete data
 N   = floor(Tf/Ts); % Number of samples for estimation
 
 % Generate sampling index points randomly
 idx = cumsum(1+rand(1,N));
 idx = [1 round(Ns*idx/max(idx))];
 
 %Generate the input (irregular sampling)
 u  = randn(1,N); 
 uc = zeros(1,Ns);
 for i=1:N,
  uc(idx(i):idx(i+1)-1) = u(i);  % Upsample the input to correct the time base
 end
 uc(end) = u(end);
end

x  = zeros(2,1);      % Initial state set to zero
yc = zeros(1,Ns);     
Qs = sqrtm(Q);

% Simulate  
for k=1:Ns,
 yc(k) = C*x + 1e-2*sqrt(K(1))*randn(1)/T;
 x     = A*x + B*uc(k) + Qs*randn(2,1);
end

% Sample the data using integration (average) the data
% IMPORTANT - make sure Z.t is set so that continuous data can be recognized

yavg = zeros(1,N);
% Assume integration as sampling step
del = min(diff(idx))-1;
for i=1:N,
 yavg(:,i) = mean(yc(idx(i):idx(i)+del),2);
end

% Build up data structure from simulation results
Z.y = yavg;
Z.u = u(:)';
Z.t = ts(idx);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify model structure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.A    = 2;
M.type = 'ss';
M.op   = 's';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedure runs
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.dsp    = dsp;
OPT.miter  = 200;
OPT.alg    = 'gn';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = est(Z,M,OPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display results 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,
 Mt.disp.colour='b';
 Mt.disp.legend = 'True Response';
 Mt.w = G.w; 
 showbode(Mt,G);
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
