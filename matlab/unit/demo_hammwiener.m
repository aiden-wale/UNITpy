%  Running this demos PEM estimatin of Hammerstein-Weiner Model Structure

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

T   = 1;        % Sampling Period in seconds
N   = 1000;     % Duration of Data Recird
var = 1e-2;     % Measurement noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Linear System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

den     = real(poly([-1,-0.1,-0.01-j,-0.01+j,-0.0001-2*j,-0.0001+2*j]));
num     = 10*den(length(den));
[bq,aq] = tfdata(c2d(tf(num,den),T,'zoh'),'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Nonlinear System Parameterisation                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_type = 'deadzone';
up_out   = 0.5; 
low_out  = -0.6;    % `Output non-linearity

in_type = 'saturation';
up_in = 0.21; 
low_in = -0.41;     % Input non-linearity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t     = 0:1:N-1; 
Z.u   = sign(sin(5*pi*t/N));
z     = sat(Z.u,low_in,up_in,1);
x     = filter(bq,aq,z); 
noise = sqrt(var)*randn(size(x));
Z.y   = dzone(x(:),low_out,up_out)+noise(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mq.in.type   = in_type;  
Mq.out.type  = out_type;  
Mq.A         = length(aq)-1; 
Mq.B         = length(bq)-2;
Mq.type      = 'oe';
Mq.delay     = 1;
Mq.T         = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.miter = 1000;  %Set the maximum number of iterations
Gq = est(Z,Mq,OPT);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp,
 Gt.A = aq;
 Gt.B = bq;
 Gt.in.type = in_type;
 Gt.in.lower = low_in;
 Gt.in.upper = up_in;
 Gt.out.type = out_type;
 Gt.out.lower = low_out;
 Gt.out.upper = up_out;
 Gt.disp.colour='b';
 Gt.disp.legend = 'True Response';
 Gq.disp.error = 0;
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









