%  Running this demos MISO linear system estimation
clear; close all;
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end
OPT.dsp = dsp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T   = 1;       % Sampling Period
N   = 1000;    % Number of Samples
var = 1e-3;    % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True Linear System
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den1      = real(poly([-0.1,-1,-0.3,-0.4]));  % First system
num1      = den1(length(den1));
den2      = real(poly([-0.2,-0.1,-0.5]));     % Second System
num2      = 1.1*den2(length(den2));
[bq1,aq1] = c2dm(num1,den1,T,'zoh');     % Discrete time versions
[bq2,aq2] = c2dm(num2,den2,T,'zoh');
cq        = [1,-0.2];  
dq        = [1,-0.5];           % Measurement noise colouring

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t    = 0:1:N-1; 
u1   = sign(sin(2*pi*t/N)); 
u2   = sign(sin(3*pi*t/N));
u = [u1(:),u2(:)];
%u = randn(N,2);  % Get 2 white noise input
noise = filter(cq,dq,sqrt(var)*randn(N,1));   % Measurement Noise
y = filter(bq1,aq1,u(:,1))+filter(bq2,aq2,u(:,2))+noise(:);
Z.y = y; Z.u = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.nA    = [length(aq1)-1;length(aq2)-1];
M.nB    = [length(bq1)-2;length(bq2)-2];
M.nC    = length(cq)-1;
M.nD    = length(cd)-1;
M.delay = [1;1];
M.T     = T;
M.type  = 'bj';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = est(Z,M,OPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,
	Gt.A = [aq1;[aq2 0]];
	Gt.B = [bq1;[bq2 0]];
	Gt.C = cq;
	Gt.D = dq; 
	Gt.T = T;
 Gt.w = G.w;
	Gt.type = 'bj';
	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';
	
	showbode(Gt,G);
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













