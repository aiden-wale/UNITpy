%  Running this demos estimation of MISO Hammerstein-Wiener 
%  model structure using prediction error method.

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

T   = 1;       % Sampling Period (sec)
N   = 1000;    % Number of Samples
var = 1e-4;    % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True Linear System Components
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
%  Specify True Nonlinear System Components
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in1_type = 'saturation';
in1_low = -0.8;  in1_up  = 0.9;

in2_type = 'deadzone';
in2_low = -0.5;  in2_up  = 0.6;

out_type = 'deadzone'; 
out_low = -0.05; out_up = 0.05; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u      = randn(N,2);  % Get 2 white noise input
x(:,1) = sat(u(:,1),in1_low,in1_up,1);
x(:,2) = dzone(u(:,2),in2_low,in2_up);
noise  = sqrt(var)*randn(N,1);   % Measurement Noise
z      = filter(bq1,aq1,x(:,1))+filter(bq2,aq2,x(:,2)) + noise(:);
y      = dzone(z,out_low,out_up);
Z.y = y; Z.u = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.A    = [length(aq1)-1;length(aq2)-1];
M.B    = [length(bq1)-2;length(bq2)-2];
M.delay = [1;1];
M.T     = T;
M.type  = 'oe';

M.in(1).type  = in1_type;
M.in(2).type  = in2_type;

M.out.type    = out_type;

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
	Gt.T = T;
	Gt.type = 'oe';
	Gt.in(1).type = in1_type;
	Gt.in(1).lower = in1_low;
	Gt.in(1).upper = in1_up;
	Gt.in(2).type = in2_type;
	Gt.in(2).lower = in2_low;
	Gt.in(2).upper = in2_up;
	Gt.out.type = out_type;
	Gt.out.lower = out_low;
	Gt.out.upper = out_up; 
	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';
	
	showbode(Gt,G);
	
	figure;
	utest = linspace(-1,1,1000);
	plot(utest(:),sat(utest(:),G.in(1).lower,G.in(1).upper,1),'xr'); hold on;
	plot(utest(:),sat(utest(:),Gt.in(1).lower,Gt.in(1).upper,1),'-','linewidth',2);
	legend({'Estimated','True'},'location','southeast');
	title('Input #1 nonlinearity')
	grid on;
	hold off;
	
	figure;
	utest = linspace(-1,1,1000);
	plot(utest(:),dzone(utest(:),G.in(2).lower,G.in(2).upper),'xr'); hold on;
	plot(utest(:),dzone(utest(:),Gt.in(2).lower,Gt.in(2).upper),'-','linewidth',2);
	legend({'Estimated','True'},'location','southeast');
	title('Input #2 nonlinearity')
	grid on;
	hold off;
 
	figure;
	utest = linspace(-1,1,1000);
	plot(utest(:),dzone(utest(:),G.out.lower,G.out.upper),'xr'); hold on;
	plot(utest(:),dzone(utest(:),Gt.out.lower,Gt.out.upper),'-','linewidth',2);
	legend({'Estimated','True'},'location','southeast');
	title('Output nonlinearity')
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














