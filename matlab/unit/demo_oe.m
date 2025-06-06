%  Running this demos PEM estimation of Output-Error model structure

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

T   = 1;     % Sampling period in seconds
N   = 100;   % Number of samples
var = 1e-1;  % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify true (linear) system                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den       = real(poly([-0.1,-1,-0.2,-0.3,-0.5,-0.05+j*3,-0.05-j*3]));
num       = 10*den(length(den));
[a,b,c,d] = tf2ss(num,den);
[a,b]     = c2d(a,b,T); 
[bq,aq]   = ss2tf(a,b,c,d,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t     = 0:1:N-1; 
Z.u   = sign(sin(3*pi*t/N));
noise = sqrt(var)*randn(size(Z.u));
Z.y   = filter(bq,aq,Z.u)+noise; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.A     = length(aq)-2; 
Mq.B     = Mq.A-1; 
Mq.T     = T; 
Mq.delay = 1;
Mq.type  = 'oe';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gq = est(Z,Mq,OPT); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,
	Gt.A = aq;
	Gt.B = bq;
	Gt.T = T;
    Gt.w = Gq.w;
	Gt.type='oe';
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

