%  Running this tests estimation of ARX model structure from frequency
%  domain data

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

fs  = 10;     % Sampling frequency in Hz
N   = 1000;   % Number of samples
var = 1e-9;   % Measurement noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True System
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den = real( poly([-0.1,-0.2,-0.02+j*1,-0.02-j*1,-0.01-j*0.1,-0.01+j*0.1]) );
num = den(length(den));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w     = logspace(-3,log10(pi*fs),N);
F     = polyval(num,j*w)./polyval(den,j*w);
noise = sqrt(var)*( randn(size(F(:))) + j*randn(size(F(:))) );
F     = F(:) + noise(:)./polyval(den,j*w).';
Z     = [F(:),w(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ms.nA   = length(den)-1; 
Ms.nB   = length(num)-1;     
Ms.op   = 's'; 
Ms.T    = 1/fs; 
Ms.type = 'arx';
Ms.w    = w(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gs = est(Z,Ms,OPT); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,
	Gt.A    = den;
	Gt.B    = num;
	Gt.type = 'arx';
	Gt.op   = 's';
	Gt.T    = 1/fs;
	Gt.w    = w(:);
	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';
	Z = startZ(Z);
	Z.disp.legend = 'Measured Response';

	showbode(Gt,Gs,Z);
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













