%  Running this demos signal estimation using the 
%  Kalman predictor, filter and smoother algorithms.

clear; close all;  
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randn('state',0); rand('state',0);

T = 1;       % Sampling period in seconds
N = 50;      % Number of samples
var = 1e-1;  % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Linear System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den = real(poly([-0.1,-1,-0.2,-0.3,-0.5,-0.05+j*3,-0.05-j*3]));
num = 10*den(length(den));
[Mq.ss.A,Mq.ss.B,Mq.ss.C,Mq.ss.D] =tf2ss(num,den); 
[Mq.ss.A,Mq.ss.B] = c2d(Mq.ss.A,Mq.ss.B,T); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0:1:N-1;
u = sign(sin(10*pi*t/N));
noise = sqrt(var)*randn(size(u));
y = Mq.ss.C*ltitr(Mq.ss.A,Mq.ss.B,u(:)).'+Mq.ss.D*u; 
Z.y=y(:)+noise(:); Z.u=u(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.T = T; Mq.op = 'q';  
Mq.ss.R = var;                       % Measurement Noise Variance
Mq.ss.Q = 0.001*eye(size(Mq.ss.A));  % State Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how 
%  estimation procedure runs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.dsp = dsp;
OPT.alg = 'sqroot'; 
OPT.allP = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G  = kf(Z,Mq,OPT);
Gs = ks(Z,Mq,OPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dsp, 
	h=plot([y(:)+noise(:),y(:),G.yf(:),G.yp(:),Gs.ys(:)]);
	grid
	title('Observed data, noise free and Kalman Predictor/Filter/Smoother Ouput')
	legend('Observed','Noise Free data','Filter','Predictor','Smoother')
	set(h,'Linewidth',2);
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




