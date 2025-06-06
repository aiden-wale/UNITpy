% Running this demo's estimation of continuous time 
% tranfer function OE structure from time domain data. 

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

fs = 50;      % Sampling frequency for simulation of true system
T  = 1/fs;    % Associated samplng period for simulation
var = 1e-2;   % Measurement noise variance   
regular = 1;  % Set to 0 for irregularly spaced time samples

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a true linear system
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den = conv([1,5],[1,6]);
num = den(end)*[1,1];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate some data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0:T:10-T;                   % Easy way - lsim will  
u = randn(size(t));
u = sign(cos(4*pi*t/t(end)));
y = lsim(num,den,u,t);
Z.y = y; Z.u = u; Z.t = t;

% Could also simulate while being a bit more careful about how the 
% sampling is assumed done (i.e. integrated sampling)

Mt.T = T;
Mt.nB = length(num)-1;
Mt.ny=1;
TH = [num(:);den(2:end)'];
Mt = t2m_soe(Mt,TH);
Mt.B=num; Mt.A=den;
m = samplek(Mt);

uy  = [Z.u(:) Z.y(:)];
xh  = ltitr(m.ss.A-m.ss.K*m.ss.C,[m.ss.B-m.ss.K*m.ss.D m.ss.K],uy,m.ss.X1);
yh  = xh*m.ss.C.' + Z.u(:)*m.ss.D.';
pe  = y-yh;

% Add some noise on top of either y or yh depending on 
% the sort of sampling assumptions you want to experiment with

Z.y = yh(:) + sqrt(var)*randn(size(yh(:))); 

% Check if user now wants irregular sampling to be demo'd
if (~regular)
 yy = []; uu = []; tt = [];
 for k=1:length(y)
  if rand<0.9  % Throw samples away with 30% chance
   yy = [yy,Z.y(k)]; uu = [uu,Z.u(k)]; tt = [tt,Z.t(k)];
  end
 end;
 Z.y = yy; Z.u = uu; Z.t = tt;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify model structure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.A=2; M.B=1; M.op='s'; M.type = 'oe';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = est(Z,M);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp
	Mt.disp.colour='b';
	Mt.disp.legend = 'True Response';
	Mt.w = G.w;
 Mt.op = 's';

 %bode(ss(Mt.ss.A,Mt.ss.B,Mt.ss.C,0),G.sysG);
 
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
