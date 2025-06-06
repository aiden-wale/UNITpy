%  Running this demos estimation of ARMAX model
%  using prediction error method

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

T   = 1;       % Sampling period in seconds
N   = 1000;    % Number of samples
var = 1e-4;    % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den     = real(poly([-0.1,-1,-0.2,-0.3,-0.5]));  % Cts time spec
num     = den(length(den));
[bq,aq] = c2dm(num,den,T,'zoh');                 % Discrete time version
cq      = [1,-0.2];  
dq      = aq;               % Measurement noise colouring

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u     = randn(1,N);
noise = filter(cq,dq,sqrt(var)*randn(size(u)));
y     = filter(bq,aq,u(:))+noise(:); 
Z.y   = y; 
Z.u   = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.T  = T;
Mq.A = length(aq)-1; 
Mq.B = length(bq)-2; 
Mq.C = length(cq)-1;
Mq.delay = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gq = est(Z,Mq,OPT); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dsp, 
	Gt.A = aq;
	Gt.B = bq;
	Gt.C = cq;
	Gt.type = 'armax';
	Gt.T = T;
 Gt.w = Gq.w; 
	Gt.disp.colour='b';

	Gt.disp.legend = 'True Response';
	Gq.disp.legend = [Gq.disp.legend, ' q operator'];

	showbode(Gt,Gq);
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









