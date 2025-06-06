%  Running this demos MISO Hammerstein system estimation

clear; close all;
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end
OPT.dsp = dsp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T   = 1;       % Sampling Period
N   = 1000;    % Number of Samples
var = 1e-4;    % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Nonlinear System Components
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

low = -0.5; 
up  = 0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u      = randn(N,2);  % Get 2 white noise input
x(:,1) = sat(u(:,1),low,up,1);
x(:,2) = dzone(u(:,2),low,up);
noise  = sqrt(var)*randn(N,1);   % Measurement Noise
y      = filter(bq1,aq1,x(:,1))+filter(bq2,aq2,x(:,2)) + noise(:);
Z      = [y(:),u];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.nA    = [length(aq1)-1;length(aq2)-1];
M.nB    = [length(bq1)-2;length(bq2)-2];
M.delay = [1;1];
M.T     = T;
M.type  = 'oe';
M.in(1).type  = 'hinge';
M.in(2).type  = 'deadzone';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gq = est(Z,M,OPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp, 
	Gt.A = [aq1;[aq2 0]];
	Gt.B = [bq1;[bq2 0]];
	Gt.T = T;
	Gt.type = 'oe';
	Gt.in(1).type = 'saturation';
	Gt.in(1).lower = low;
	Gt.in(1).upper = up;
	Gt.in(2).type = 'deadzone';
	Gt.in(2).lower = low;
	Gt.in(2).upper = up;
	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';
	
	showbode(Gt,Gq);
	
	figure;
	utest = linspace(-1,1,1000);
	plot(utest(:),hinge(utest(:),Gq.in(1).eta),'xr'); hold on;
	plot(utest(:),sat(utest(:),Gt.in(1).lower,Gt.in(1).upper,1),'-','linewidth',2);
	legend({'Estimated','True'},'location','southeast');
	title('Input #1 nonlinearity')
	grid on;
	hold off;
	
	figure;
	utest = linspace(-1,1,1000);
	plot(utest(:),dzone(utest(:),Gq.in(2).lower,Gq.in(2).upper),'xr'); hold on;
	plot(utest(:),dzone(utest(:),Gt.in(2).lower,Gt.in(2).upper),'-','linewidth',2);
	legend({'Estimated','True'},'location','southeast');
	title('Input #2 nonlinearity')
	grid on;
	hold off;
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














