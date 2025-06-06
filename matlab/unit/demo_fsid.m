%  Running this demos Frequency Domain Subspace-based estimation of state
%  space model.

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

fs =  2;      % Sampling Frequency
N  =  2000;   % Number of samples
var = 1e-4;  % White measurement noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True Linear System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den     = real( poly([-0.02+j*1,-0.02-j*1,-0.01-j*0.1,-0.01+j*0.1]) );
%den     = real( poly([-0.1,-0.2,-0.02+j*1,-0.02-j*1,-0.01-j*0.1,-0.01+j*0.1]) );
num     = 10*den(length(den));
sys     = tf(num,den); 
sysd    = c2d(sys,1/fs,'zoh');
[bq,aq] = tfdata(sysd,'v');
[bd,ad] = s2delta(num,den,1/fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w          = logspace(-3,log10(pi*fs),N);
F          = polyval(num,j*w)./polyval(den,j*w);
noise      = sqrt(var)*( randn(size(F(:))) + j*randn(size(F(:))) );
F          = F(:) + noise(:); 
Z.y(1,1,:) = F(:);  
Z.w        = w(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ms.w = w; 
Mq.w = w; 

Ms.nx = length(den)-1; 
Ms.op = 's'; 
Ms.T = 0; 

Mq.nx = length(den)-1; 
Mq.op = 'q'; 
Mq.T = 1/fs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gs = fsid(Z,Ms);
Gq = fsid(Z,Mq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,
	%Gt = polyval(num,j*w)./polyval(den,j*w);
	Gt.G = polyval(num,j*w(:))./polyval(den,j*w(:));
	Gt.w = Gq.w;
	Data.G = F(:)+noise(:);
	Data.w = Gt.w;

	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';
	Gq.disp.legend = [Gq.disp.legend, ' q operator'];
	Gs.disp.legend = [Gs.disp.legend, ' s operator'];
	Data.disp.legend = 'Noise corrupted data';
	Gt.disp.aux='magonly';

	showbode(Gt,Gs,Gq,Data);
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













