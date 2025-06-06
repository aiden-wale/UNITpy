%  Running this demos non-parametric estimation algorithms
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

T   = 1;      % Sampling Period (sec)
N   = 20000;   % Number of samples
var = 1e-4;   % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True Linear System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den = real(poly([-0.1,-1,-0.2,-0.3,-0.5]));  % Cts time spec
num = den(length(den));
[bq,aq] = c2dm(num,den,T,'zoh');             % Discrete time version

cq = [1,-0.2];  dq = [1,-0.5];               % Measurement noise colouring

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = randn(1,N);
k = 0:1:N-1; u = sign(cos(2*pi*k*6/N)); 
u = [u(1:floor(N/2)),randn(1,floor(N/2)+1)]; 
noise = filter(cq,dq,sqrt(var)*randn(size(u)));
y = filter(bq,aq,u); 
Z.y = y(:)+noise(:); Z.u = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.type = 'nonpar'; 
M.T = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gbltk = est(Z,M,OPT); 

OPT.alg = 'etfe';
Getfe = est(Z,M,OPT); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
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
	
	showbode(Gt,Gbltk,Getfe);
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









