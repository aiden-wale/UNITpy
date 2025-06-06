%  Running this demo's ML estimation of MIMO systems from frequency
%  domain data using a subspace method, gradient based search, and the 
%  the EM algorithm.

clear; close all;
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 200;  % Number of data points
T=1/500;  % Sample time

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a true linear system
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%System
% n=6; m=2; p=2;
% sys=drss(n,p,m); 
% [a,`b,c,d]=ssdata(sys); 
% a=0.999*a;
% sys=ss(a,b,c,d,-1);

% n=6; p=1; m=1;
% den = real( poly([-0.1,-0.2,-0.02+j*1,-0.02-j*1,-0.01-j*0.1,-0.01+j*0.1]) );
% num = 10*den(length(den));
% [a,b,c,d] =tf2ss(num,den); 
% sysc=ss(a,b,c,d); 
% sys=c2d(sysc,1,'zoh');
% [a,b,c,d]=ssdata(sys);

n=2; m=2; p=2;
eta=0.001;
omegad=logspace(1,2,n);
omegan=omegad+0.1*omegad;
num=1; den=1;
for i=1:n,
 num=conv(num,[1 2*eta*omegan(i) omegan(i)^2]);
 den=conv(den,[1 2*eta*omegad(i) omegad(i)^2]);
end
if m>1 || p>1,
 nums = num; clear num;
 dens = den; clear den;
 for i=1:m,
  for j=1:p,
   num{j,i,:} = nums + 0.1*randn(size(nums)); 
   den{j,i,:} = dens;
  end
 end  
end

w=logspace(-2,log10(pi/T),N); 
sysc=tf(num,den);
sysd=c2d(sysc,T,'zoh');
sys=ss(sysd);
[a,b,c,d]=ssdata(sys);
n=size(a,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=zeros(n,m,N); v=zeros(n,m,N); y=zeros(p,m,N); e=zeros(p,m,N); ejw=exp(sqrt(-1)*w*T);
for k=1:N,
 x(:,:,k)=(eye(n)*ejw(k)-a)\b;  
 e(:,:,k)=0.1*(randn(p,m)+sqrt(-1)*randn(p,m));
 y(:,:,k)=c*x(:,:,k) + d + e(:,:,k); 
   
 ynf(:,:,k)=c*((eye(n)*ejw(k)-a)\b) + d; %no noise 
end

z.y = y; 
z.w = w(:); 
z.T = T;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify model structure
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q=1*eye(n); R=1e-1*eye(p);   %Initial guess at covariances

mm.A    = n; 
mm.op   = 'q'; 
mm.T    =T; 
mm.type ='ss';
mm.w    =w(:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Optional parts about how the estimation procedure runs
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oss.alg='sid';              % Options for subspace estimation
oss.lag=round((N-10)/2);
oss.dsp=dsp;

ogn.dsp   = dsp;            % Options for GN-search estimation
ogn.par   = 'ddlc'; 
ogn.cost  = 'det'; 
ogn.op    = 'q';
ogn.dir   = 'trust';
ogn.miter = 100; 
ogn.ngt   = 0;

oem.dsp     = dsp;          % Options for EM-alg estimation
oem.miter   = 200; 
oem.optit   = 100;
oem.alg     = 'em';
oem.stoptol = 1e-5; 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gss = est(z,mm,oss);     % Subspace-based estimate
ggn = est(z,gss,ogn);    % ML via GN-search estimate starting at sid estimate

gss.ss.Q = Q;            % Reset Q and R matrices to 
gss.ss.R = R;            % Initial values
gem = est(z,gss,oem);    % ML EM search starting at sid estimate

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Display the results
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dsp, 
 gt.ss.A=sys.A; 
 gt.ss.B=sys.B; 
 gt.ss.C=sys.C; 
 gt.ss.D=sys.D; 
 gt.ss.K=[];
 gt.T=T; 
 gt.type='ss'; 
 gt.op='q'; 
 gt.w=w(:);
 gt.disp.colour = 'b';
 gt.disp.legend = 'True Response';
 gt.disp.aux    = 'magonly';
 gt = startM(z,gt);
 data.G=y; 
 data.w=w(:); 
 data.disp.legend = 'Data';
 
 showbode(gt,data,gss,gem,ggn);
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
