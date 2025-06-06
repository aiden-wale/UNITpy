%  Running this demos estimation of ARX model structure 
%  using recursive least squares algorithm

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

T = 0.1;    % Sampling Period (sec)
N = 100;    % Number of Samples
var = 1e-3; % Measurement Noise Variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify true linear System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den = real(poly([-0.1,-1,-0.2]));
num = den(length(den));
[bq,aq] = c2dm(1,den,T,'zoh');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = randn(1,N);
noise = sqrt(var)*randn(size(u));
y = filter(bq,aq,u); 
noise = filter(100*sum(aq),aq,noise(:));
Z = [y(:)+noise(:),u(:)]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.w = logspace(log10(0.001/T),log10(pi/T),1000); Md.w = Mq.w;
Mq.A = length(aq)-1; Mq.B = Mq.A; Mq.T = T; Mq.lnorm = 1; Mq.type = 'arx';
Md = Mq; Mq.op='q'; Md.op='d';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how 
%  estimation procedure runs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.dsp = dsp;
OPT.n=0; 
OPT.alg.type = 'rls'; 
OPT.alg.lambda = 0.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp, disp('Computing Shift Operator Estimate...'); end
Gq = est(Z,Mq,OPT);
if dsp, disp('Computing Delta Operator Estimate...'); end
Gd = est(Z,Md,OPT); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp, 
	ww = exp(j*Mq.w*T);
	Gt = polyval(bq,ww)./polyval(aq,ww);
	
	subplot(211)
	h=semilogx(Mq.w,log10(abs([Gt(:),Gq.G(:),Gd.G(:)])));
	set(h,'Linewidth',2);
	legend('G (true)','G(q) (est)','G(d) (est)')
	grid
	title('True and Estimated Frequency Responses')

	subplot(212)
	h=plot(Gq.th_hist);
	set(h,'Linewidth',2);
	grid
	title('Time Evolution of Shift Operator Estimates')
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





