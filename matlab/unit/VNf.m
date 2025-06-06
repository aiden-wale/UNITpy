%  Function to calculate quadratic cost
%
%  V_N(theta) = 1/N\sum_{k=1}^N[F(w_k) - B(p)/A(p)cos(2*pi*w_k)]^2
%
%  Associated with output error model for observed frequency domain data:
%
%  F(w_s) = B(p)/A(p)cos(2*pi*w_s)
%
%  The operator p can be the Z tranform variable with z=e^(j*w_s*T) (with T
%  being the sampling period in seconds), the Euler differencing (delta)
%  operator d = (q-1)/T with d = (e^(j*w_s*T)-1)/T or the Laplace Transform
%  variable s with s=j*w.
%
%  This function is not meant to be directly called by users - instead it
%  is an auxiliary function used by foe.m for frequency domain system
%  identification.
%
%  Usage is
%
%  cost = VNf(Z,theta,OPT.n)
%
%  Z        = observed frequency response data [F(:),w(:)] where
%             plot(w,abs(F)) should plot the measured  frequency
%             response.   Units for w are real *not* normalised freq.
%  theta    = [b(d),a(d)]: specification of point to calculate cost at.
%  M        = Data structure which defines the model structure which
%             is to be estimated from the data as follows:
%  OPT      = Data structure which defines options for an estimation
%             algorithm.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.


% Copyright (C) Brett Ninness.

function [cost,pe,grad,phi] = VNf(Z,theta,OPT,M,div)

if nargin < 5 div = 0; end;

% Extract out relevant vectors from input data
[F,w,ny,nu,Ny] = Z2data(Z); F=squeeze(F); F=F(:); wmax = 1;%max(w);

%Make input equal to the identity for each w(k) if not supplied by user
if ~isfield(Z,'u'),
 Z.u = zeros(1,1,N);
 for k=1:N, Z.u(:,:,k) = 1; end
elseif isempty(Z.u),
 Z.u = zeros(1,1,N);
 for k=1:N, Z.u(:,:,k) = 1; end
end
U = squeeze(Z.u); U = U(:);

if ~exist('OPT')
 OPT.W     = ones(size(F));
else
  if ~isfield(OPT,'W') OPT.W=ones(size(F));      end;
end;

%  Determine frequency domain argument according to time domain operator
if (M.op=='q'), 
 z = exp(j*w*M.T);
elseif (M.op == 's')
 z = j*w;
elseif (M.op == 'd') 
 z = (exp(j*w*M.T)-ones(size(w)))/M.T;
end

%  Extract numerator and denominator from theta parameter vector.
b = theta(1:M.nB+1); 
a = theta(M.nB+2:end); 
m = length(b);
n = length(a)+1;
mn = max(m,n);

PHI = zeros(length(w),mn);
X   = eye(mn,mn);
%  Then use it to generate regressors
for k=1:mn 
 PHI(:,k) = polyval([X(mn-k+1,:)],z);
end;
PHI = fliplr(PHI);  % order L->R from lowest order to highest order poly

%Extract weighting vector
weights = sqrt(OPT.W(:));

% Calculate cost and return it
num  = PHI(:,end-m+1:end)*b; 
den  = PHI(:,end-n+1:end)*[1;a]; 
fhat = (num./den).*U(:);
pe   = F(:)-fhat(:); 
pe   = weights.*pe;
cost = 0.5*real(pe'*pe)/length(pe);

% Calculate gradient and Hessian if requested
if div
 PSI = zeros(length(w),length(theta));
 for k=1:m
  nb = [zeros(1,k-1),1,zeros(1,m-k)];
  PSI(:,k) = weights.*((PHI(:,end-m+1:end)*nb(:))./den).*U(:);
 end;

 for k=2:n
  nb = [zeros(1,k-1),1,zeros(1,n-k)];
  PSI(:,k+m-1) = -weights.*(PHI(:,end-n+1:end)*nb(:))./(den).*fhat;
 end;

 npar=size(PSI,2);
 R=triu(qr([[real(PSI);imag(PSI)], [real(pe);imag(pe)]]));
 phi=R(1:npar,1:npar)/sqrt(Ny); 
 pe=-R(1:npar,end)/sqrt(Ny);
 grad=phi'*pe;
end;




