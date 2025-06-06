%  Running this demos signal estimation using Sequential 
%  Important Re-Sampling (SIR) - also known as particle 
%  filtering.  This is done for a linear Gaussian system 
%  the approximate SIR filter can be compared to the exact
%  Kalman filter answer.

clear; close all;  
global dm;  if isempty(dm),  clear global dm;  dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 1;       % Sampling Period
N = 50;      % Number of samples
var = 1e-2;  % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Linear System - randomly drawn                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randflag = 0;  % Set to 1 if you want a randomly chosen system

if randflag
 nx=10; nu=1; ny=1; sysd=drss(nx,ny,nu);
 [M.ss.A,M.ss.B,M.ss.C,M.ss.D] = ssdata(sysd); 
else
 den = real(poly([-0.1,-1,-0.2,-0.3,-0.5,-0.05+j*3,-0.05-j*3]));
 num = 10*den(length(den));
 [M.ss.A,M.ss.B,M.ss.C,M.ss.D] =tf2ss(num,den); 
 [M.ss.A,M.ss.B] = c2d(M.ss.A,M.ss.B,T); 
 nx = size(M.ss.A,1); nu = 1; ny = 1;
end;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if randflag
 u = randn(nu,N); t = 0:1:N-1;
 y = lsim(sysd,u);
 noise = sqrt(var)*randn(N,ny);
 Z.y = y(:)+noise(:); Z.u = u;
else  
 t = 0:1:N-1;
 u = sign(sin(10*pi*t/N)); 
 y = M.ss.C*ltitr(M.ss.A,M.ss.B,u(:)).'+M.ss.D*u; 
 noise = sqrt(var)*randn(N,ny);
 Z.y=y(:)+noise(:); Z.u=u(:);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.T = T; M.op = 'q';  
M.ss.R = var;             % Measurement Noise Variance
M.ss.Q = 0.001*eye(nx);   % State Noise variance
M.ss.X0= zeros(nx,1);      % Prior mean on initial state 
M.ss.P0= eye(nx);        % Covariance on initial state

M.model=@ssmod;           % Handle to function that specifies
                          % model for particle filter (SIR) alg.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how 
%  estimation procedure runs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.alg  = 'sqroot';    % Use square root forms of algorithms
OPT.allP = 0;           % Don't store and hand back state covariances
OPT.pnum = 500;         % Number of particles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = sir(Z,M,OPT);   % Approximate SIR-based predictor and filter
G = kf(Z,M,OPT);    % Exact Kalman predictor and filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp
 h=plot([Z.y(:),y(:),G.yf(:),g.yf(:)]);
 grid
 title('Observed data versus Kalman Predictor, and Particle predictor')
 legend('Observed','Noise Free data','Kalman Filter','Particle Filter')
 set(h,'Linewidth',2);   
end;

if dm
 disp('  ')
 disp('---------------------------------------------------------------------')
 disp('  ')
 disp('You now have access to the MATLAB workspace so that you may examine')
 disp('the results of this simulation.  To return to the demos, type "dbcont"')
 disp(' ')
 keyboard; 
end;




