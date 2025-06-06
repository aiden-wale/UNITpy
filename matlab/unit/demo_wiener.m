%  Running this demos estimation of Weiner Model
%  via prediction error method computed by gradient
%  based search algorithm

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

T   = 1;        % Sampling Period (seconds)
N   = 100;      % Duration of Data Record (# of samples)
var = 1e-3;        % Measurement noise innovations variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True Linear System  Component                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den     = real(poly([-1,-0.5,-0.01-0.8*j,-0.01+0.8*j,-0.0001-1.0*j,-0.0001+1.0*j]));
num     = 0.4*den(length(den));
[bq,aq] = tfdata(c2d(tf(num,den),T,'zoh'),'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True Nonlinear System Component                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = 'deadzone';
up = 0.21; low = -0.41;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z.u   = randn(1,N);   % Type of input
t     = 0:1:N-1; 
Z.u   = sign(sin(3*pi*t/N));
x     = filter(bq,aq,Z.u(:));
noise = sqrt(var)*randn(size(Z.u));
Z.y   = dzone(x(:),low,up)+noise(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.out.type  = type;  % Output Nonlinearity Type
Mq.A = length(aq)-1; 
Mq.B = length(bq)-2;
Mq.delay = 1;
Mq.T = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gq = est(Z,Mq,OPT);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,
	Gt.A = aq;
	Gt.B = bq;
	Gt.type = 'oe';
	Gt.T = T;
	Gt.out.type = 'deadzone';
	Gt.out.lower = low;
	Gt.out.upper = up;
	Gt.disp.colour='b';
	Gt.disp.legend = 'True Response';
	Gq.disp.error = 0;
		
	showbode(Gt,Gq);
	
	figure;
	utest = linspace(-1,1,1000);
	plot(utest(:),dzone(utest(:),Gq.out.lower,Gq.out.upper),'xr'); hold on;
	plot(utest(:),dzone(utest(:),Gt.out.lower,Gt.out.upper),'-','linewidth',2);
	legend({'Estimated','True'},'location','southeast');
	title('Output nonlinearity')
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









