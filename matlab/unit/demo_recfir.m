%  Running this demos estimation of FIR model structure
%  using recursive least squares algoritm

clear; close all;  
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 0.1;    % Sampling period in seconds
N=1000;     % Number of data sample to generate
var = 1e-4; % Measurement noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify true linear System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bq = impulse(c2d(tf(1,real(poly([-1 -3+j -3-j]))),0.1,'zoh'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u     = randn(1,N);
noise = sqrt(var)*randn(size(u));
y     = filter(bq,1,u)+noise; 
Z.y   = y; 
Z.u   = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.w    = logspace(log10(0.001/T),log10(pi/T),1000); 
Mq.T    = T; 
Mq.nB   = length(bq)-1;
Mq.type = 'fir';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how 
%  estimation procedure runs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.dsp = dsp;
OPT.n=0; 
OPT.alg.type = 'rls';

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
	ww = exp(j*Mq.w*T);
	Gt = polyval(bq,ww)./polyval(1,ww);
	subplot(211)
	h=semilogx(Mq.w,20*log10(abs([Gt(:),Gq.G(:)])));
	set(h,'Linewidth',2);
	title('True and Estimated Frequency Responses')
	grid
	legend('Gt','Gq');
	subplot(212)
	h=plot(Gq.th_hist)';
	set(h,'Linewidth',2);
	grid
	title('Time evolution of Parameter Estimates')
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

clear OPT;






