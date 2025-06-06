%  Running this demos least squares estimation of an 
%  orthonormal basis model structure

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

T   = 0.1;    % Sampling period in seconds
N   = 1000;   % Number of data sample to generate
var = 1e-31;  % Measurement noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify true (linear) System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den = real(poly([-0.1,-1,-0.2,-0.3+j,-0.3-j]));
%den = real(poly([-0.1]));
num = den(length(den));
[bq,aq] = c2dm(num,den,T,'zoh');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u     = randn(1,N);
noise = sqrt(var)*randn(size(u));
y     = filter(bq,aq,u); 
Z     = [y(:)+noise(:),u(:)]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.w     = logspace(log10(0.001/T),log10(pi/T),1000); 
Mq.A     = aq; 
Mq.poles = cplxpair(roots(aq)); 
Mq.T     = T; 
Mq.type  = 'fir';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how 
%  estimation procedure runs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.n   = 0;    
OPT.dsp = dsp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gq = onid(Z,Mq,OPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dsp,
	ww = exp(j*Mq.w(:)*T);
	Gt.G = polyval(bq,ww(:))./polyval(aq,ww(:));

	Gt.w = Gq.w; Gt.T=Gq.T;  Ghat1.w = Gq.w; Ghat1.T=Gq.T;
	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';
	Gq.disp.legend = [Gq.disp.legend, ' q operator'];
	Gt.disp.aux='magonly';

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








