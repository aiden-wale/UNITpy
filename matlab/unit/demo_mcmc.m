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

T = 1;                  % Sampling Period
N = 20;                 % Number of samples
var = 1e-2;             % White Measurement Noise variance
ueps = sqrt(3*var);     % Equivalent Uniform density bounds
dens = 'uniform';       % Could be `gaussian' or 'uniform';
w = logspace(-3,pi,500);
ww = exp(j*w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify True Linear System                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mtype = 1;
switch mtype
 case 1  % Simple first order system
  delay = 0;
  aq = [1,-0.8];
  bq = [zeros(1,delay),sum(aq)];
  cq = [];
  dq = [];
 case 2  % Simple second order system
  delay = 0;
  aq = poly([0.8,0.5]);
  bq = [zeros(1,delay),sum(aq)/4 3*sum(aq)/4];
  cq = [];
  dq = [];
 case 3  % Resonant 2nd order system
  delay = 0;
  xi = [0.95*exp(j*pi/3),0.95*exp(-j*pi/3)]; % z-domain poles
  aq = real(poly(xi)); 
  bq = [1,-0.5];
  bq = bq*sum(aq)/sum(bq);
  cq = [];
  dq = [];        
 case 4  % Resonant 4th order system
  delay = 0;
  bq = poly([-8.0722,-0.8672,0.0948]);
  aq = real(poly([0.75*exp(j*pi/3),0.75*exp(-j*pi/3),0.95*exp(j*pi/12),0.95*exp(-j*pi/12)])); 
  aq = real(poly([0.99*exp(j*pi/3),0.99*exp(-j*pi/3),0.99*exp(j*pi/12),0.99*exp(-j*pi/12)]));   
  bq = bq*sum(aq)/sum(bq);
  cq = [];
  dq = [];        
 case 5  % Non-resonant 4th order system 
  delay = 0;
  bq = [0,0,1];
  aq = real(poly([0.15*exp(j*pi/12),0.15*exp(-j*pi/12),0.5,0.8])); 
  bq = bq*sum(aq)/sum(bq);
  cq = [];
  dq = [];        
 case 6 % Resonant 4th order system used in paper
  delay = 0;
  xi = [0.75*exp(j*pi/3),0.75*exp(-j*pi/3),0.95*exp(j*pi/12),0.95*exp(-j*pi/12)]; 
  den = real(poly(log(xi))); 
  num = den(length(den));
  [bq,aq] = c2dm(num,den,1,'zoh');
  cq =[]; 
  dq = [];
 case 7  % Åström system
  delay = 0;
  aq = [1,-1.5,0.7];  
  bq = [1,0.5];  bq = bq*sum(aq)/sum(bq);
  cq = [];
  dq = [];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = randn(1,N);
u = ones(1,N); u(floor(N/2):end) = zeros(1,length(u(floor(N/2):end)));
if strcmp(lower(dens),'gaussian')
 noise = sqrt(var)*randn(size(u));
elseif strcmp(lower(dens),'uniform')
 noise = 2*ueps*(rand(size(u))-0.5*ones(size(u)));
end;
y = filter(bq,aq,u);  
Z.y = y(:)+noise(:);  Z.u=u;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Model Structure                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M.w = logspace(-3,pi,1000);  M.T=T; M.op='q'; 
M.A = length(aq)-1; M.B = length(bq)-1; M.delay=delay;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how 
%  estimation procedure runs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT.dsp  = dsp;    % Display the results
OPT.Mmax = 1e4;    % Run Mmax Metropolis-Hastings iterations
OPT.dens = dens;   % Tell the MH method the true density
OPT.n=0;           % Throw away OPT.n initial data samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate a model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gest = est(Z,M,OPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Compute Posterior Distributions using MCMC methods and plot them
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ginit = gest;                % Initialise Markov Chain iterations at PEM estimate

OPT.var = 10*var;
OPT.sampler = 'slice';       % Selects slice sampler
                             %OPT.sampler = 'metropolis'; % Selects
                             %Metropolis-Hastings sampler



G = postdist(Z,ginit,OPT);   % Run MCMC algorithm to get samples from posterior

if dsp
 Gt.A = aq;
 Gt.B = bq;
 Gt.T = T;
 Gt.w = gest.w;
 Gt.disp.legend = 'True Response';
 
 shownyq(Gt,gest,G); 
 
 figure(2)
 showdist(G);                  % Plot posteriors as smoothed versions of
 
 if (mtype==1)  %  For the simplest case we will profile true values vs posteriors
  figure(2)
  hold on
   plot([aq(2),aq(2)],[0,max(G.pa.p)],'color','red','linewidth',2);
   text(1.08*aq(2),0.95*max(G.pa.p),'True Value','color','red','FontSize',16)   
  hold off
  figure(3)
  hold on
   plot([bq(1),bq(1)],[0,max(G.pb.p)],'color','red','linewidth',2);
   text(1.08*aq(2),0.95*max(G.pb.p),'True Value','color','red','FontSize',16)     
  hold off
 end;
 
 
 [h.p,h.x]=hist(G.varlog,113); % Get sample histograms of measurement
                               % noise variance
 hest = kde(h);                % Smooth it via kernel density estimation
 
 figure               
 plot(hest.x,hest.p,'linewidth',2);
 title('Posterior distribution for noise variance')
 hold on;
 plot([var,var],[0,max(hest.p)],'color','red','linewidth',2);
 hold off;
 text(1.08*var,0.95*max(hest.p),'True Value','color','red','FontSize',16)
end;                         


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
