%  Running this tests PEM estimation of Box-Jenkins model structure

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
T   = 1;     % Sampling period in seconds
N   = 500;   % Number of samples
var = 1e-2;  % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den     = real(poly([-0.1,-1,-0.001-0.2j,-0.001+0.2j,-0.5]));  % Cts time spec
num     = den(length(den));
[bq,aq] = c2dm(num,den,T,'zoh');  % Discrete time version
cq      = [1,-0.2];  
dq      = [1,-0.5];               % Measurement noise colouring

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t     = 0:1:N-1; 
Z.u   = sign(sin(5*pi*t/N));
%Z.u     = randn(1,N);
noise   = filter(cq,dq,sqrt(var)*randn(size(Z.u)));
Z.y     = filter(bq,aq,Z.u)+noise; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.A = length(aq)-1; 
Mq.B = length(bq)-2;  
Mq.C = length(cq)-1; 
Mq.D = length(dq)-1; 
Mq.T = T;
Mq.delay = 1;
Mq.type  = 'bj';

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
	Gt.D = dq;
	Gt.T = T;
	Gt.type = 'bj';
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









