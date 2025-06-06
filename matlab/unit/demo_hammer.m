%  Running this demos PEM estimation of Hammerstein model structure

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

T   = 1;        % Sampling period in seconds
N   = 1000;     % Duration of record
var = 1e-2;     % Measurement noise variance            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True Linear System Component                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den     = real(poly([-1,-0.5,-0.01-j,-0.01+j,-0.0001-2*j,-0.0001+2*j]));
num     = den(length(den));
[bq,aq] = tfdata(c2d(tf(num,den),T,'zoh'),'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Ture Nonlinear System Component                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = 'saturation';
low = -0.6; up = 0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t     = 0:1:N-1; 
Z.u   = sign(sin(5*pi*t/N));
x     = sat(Z.u,low,up,1);
noise = sqrt(var)*randn(size(Z.u));
Z.y   = filter(bq,aq,x(:))+noise(:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.in.type  = type;
Mq.nA       = length(aq)-1; 
Mq.nB       = Mq.nA-1;
Mq.delay    = 1;
Mq.T        = T;

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
	Gt.type = 'oe';
	Gt.in(1).type='saturation';
	Gt.in(1).lower = low;
	Gt.in(1).upper = up;
	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';
	Gt.w = Gq.w;
 
	showbode(Gt,Gq);
	
	figure;
	utest = linspace(-1,1,1000);
	plot(utest(:),sat(utest(:),Gq.in.lower,Gq.in.upper,1),'xr'); hold on;
	plot(utest(:),sat(utest(:),Gt.in.lower,Gt.in.upper,1),'-','linewidth',2);
	legend({'Estimated','True'},'location','southeast');
	title('Input nonlinearity')
	grid on;
	hold off;
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








