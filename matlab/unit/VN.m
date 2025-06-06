%  Function to calculate quadratic cost
%
%  V_N(theta) = 1/N\sum_{t=1}^N[y_t - {(1-D/C)y_t + DB/CA u_t}]^2
%
%  associated with the Box-Jenkins  model:
%
%  y_t = B(p)/A(p)u_{t-delay} + C(p)/D(p)e_t
%
%  where e_t is white noise, p can be the backward shift operator q^{-1} or
%  the Euler differencing (delta) operator d = (q-1)/T (with T being
%  the sampling period).
%
%  Usage is
%
%  [cost,pe,grad,hes,PSI] = VN(Z,theta,OPT,M,div);
%
%  where
%
%   Z         = Input output data in the form Z = [y,u] where y is a column
%               vector of output measurements and u is a matrix whose
%               columns are the input measurements - this means that MISO
%               models are catered for, for MIMO the user should conduct
%               multiple MISO estimation runs (one for each output).
%  theta      = [b,a,c,d] = specification of point to calculate grad/hess
%  M          = Data structure which defines the model structure which
%               is to be estimated from the data as follows:
%    M.A,M.B  = Initial guess for input-output dynamics.
%    M.C,M.D  = Initial guess for measurement noise model.  If not
%               specified, the default is M.C/M.D=1;
%    M.delay  = Number of delays to include (see above model);
%    M.op     = set to 'q' for shift and 'd' for delta.  Default = 'q'.
%    M.T      = sampling period in s. (Ignored for q case) Default = 1;
%  OPT        = Data structure which defines options for the estimation
%               algorithm as follows:
%    OPT.dsp  = optional, set to 'trace' for verbose output.
%    OPT.n    = number of starting data points to discard to get
%               rid of initial condition effects.
%    OPT.M    = Model structure about which the parameter estimate is
%               regularised with weight OPT.delta.
%    OPT.delta= Regularisation weight (see above).
%    OPT.miter= Maximum number of updates of estimate from initial guess.
%    OPT.tol  = Expected improvement must by > OPT.tolx100% to continue.
%    OPT.lmax = Maximum number of times search distance will be shortened
%               by bisection.
%    OPT.step = Number of samples ahead to use in prediction error.
%  div        = flag variable, that if =1 causes gradients and hessians
%               to be calculated, but not otherwise
%
%    cost     = value of quadratic cost V_N(M).
%    pe       = prediction error sequence.
%    grad     = gradient of V_N(M).
%    hes      = estimate of Hessian of V_N(M).
%    PSI      = matrix with columns being prediction error gradients.
%
%   written by Brett Ninness  School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness

function [cost,pe,grad,phi] = VN(Z,theta,OPT,M,div)

if (nargin<5)  div = 0; end;  % Default is don't compute gradients

% Extract input and output from data matrix
[y,u,ny,nu,Ny] = Z2data(Z);

%Make sure the number of inputs and outputs comes from the model, not the
%data
nu = M.nu;
ny = M.ny;

% Include delays specified in model structure on inputs
for r=1:nu, 
 u(:,r) = [zeros(M.delay(r),1);u(1:Ny-M.delay(r),r)]; 
end

% Convert from stacked parameter vector to model structure form
Mnew = theta2m(theta,M,1);

%Extract polynomial coefficients and dimensions (which depend on type)
a = Mnew.A; 
b = Mnew.B; 
c = Mnew.C; 
d = Mnew.D;

% Apply any specified input non-linearity
if nu >0, [x,z]   = u2x(u,Mnew); end          % x=X(u_{t-k}), z=dX(u_{t-k})/d(eta).
Mout    = Mnew;
Mout.in = Mout.out;  % Set up defininitions of output non-linearity

% Calculate prediction errors according to model structure type
switch M.type,

 case {'arx','narx'}
  OPT.dsp   = 0; 
  OPT.fast  = 1;
  g         = barx(Z,Mnew,OPT); 
  [yhat,zz] = u2x(g.phi*g.th,Mout); 
  pe        = y(:)-yhat(:);

 case {'arma'}
  pex = zeros(Ny,1);

  % Now figure out appropriate filter for k=OPT.step ahead prediction.
  % First step is to break H(q) = Hbar_k(q) + H_tilde_k(q) where Hbar_k is FIR of order k=OPT.step.
  ex    = zeros(1,length(a)-length(c));
  [q,r] = deconv([c,ex],a);
  Hbar  = q;
  nc = length(c);
  for k=2:OPT.step,
   [q,r] = deconv([r(2:nc),0],a);
   Hbar  = [Hbar,q]; % Get first OPT.step terms in impulse response of H
  end;
  Hden = conv(c,[1,zeros(1,OPT.step-1)]);
  Hnum = conv(a,Hbar);

  %Compute prediction error
  ygu  = y;    % Prediction error not accounting for noise model
  pe   = ufilter(Hnum,Hden,ygu,Mnew);
 
 case {'fir','nfir'}
  OPT.dsp   = 0; 
  OPT.fast  = 0;
  g         = fir(Z,Mnew,OPT); 
  [yhat,zz] = u2x(g.phi*g.th,Mout); 
  pe        = y(:)-yhat(:);
  
 case {'armax','narmax'}
  %  Now Figure out what G(p)X(u_{t-k}) is for all the columns of u
  for k=1:nu,
   pex(:,k) = ufilter(b(k,:),c,x(:,k),Mnew);
  end

  % Now figure out appropriate filter for k=OPT.step ahead prediction.
  % First step is to break H(q) = Hbar_k(q) + H_tilde_k(q) where Hbar_k is FIR of order k=OPT.step.
  ex    = zeros(1,length(a)-length(c));
  [q,r] = deconv([c,ex],a);
  Hbar  = q;
  nc = length(c);
  for k=2:OPT.step,
   [q,r] = deconv([r(2:nc),0],a);
   Hbar  = [Hbar,q]; % Get first OPT.step terms in impulse response of H
  end;
  Hden = conv(c,[1,zeros(1,OPT.step-1)]);
  Hnum = conv(a,Hbar);

  %Compute prediction error
  pe   = ufilter(Hnum,Hden,y,Mnew) - sum(pex,2);
  
 otherwise  %OE, BJ, ARMA and nonlinear versions
  %  Now Figure out what G(p)X(u_{t-k}) is for all the columns of u
  if nu>0 && ~strcmp(lower(M.type),'arma'),
   for k=1:nu,
    pex(:,k) = ufilter(b(k,:),a(k,:),x(:,k),Mnew);
   end
  else
   pex = zeros(Ny,1);
  end
  
  % Now take output non-linearity into account
  % Get total contribution X(G_1u_1 + G_2u_2 + ...) VN
  if div, [spex,zz,w] = u2x(sum(pex,2),Mout); else spex = u2x(sum(pex,2),Mout); end;

  % Now figure out appropriate filter for k=OPT.step ahead prediction.
  % First step is to break H(q) = Hbar_k(q) + H_tilde_k(q) where Hbar_k is FIR of order k=OPT.step.
  ex    = zeros(1,length(d)-length(c));
  [q,r] = deconv([c,ex],d);
  Hbar  = q;
  nc = length(c);
  for k=2:OPT.step,
   [q,r] = deconv([r(2:nc),0],d);
   Hbar  = [Hbar,q]; % Get first OPT.step terms in impulse response of H
  end;
  Hden = conv(c,[1,zeros(1,OPT.step-1)]);
  Hnum = conv(d,Hbar);

  %Compute prediction error
  ygu  = y-spex;    % Prediction error not accounting for noise model
  pe   = ufilter(Hnum,Hden,ygu,Mnew);
end

% Calculate least-squares cost.
cost = pe(OPT.n+1:length(pe))'*pe(OPT.n+1:length(pe))/length(pe(OPT.n+1:length(pe)));

% Calculate gradient and Hessian of cost if it is requested by setting flag div
if div
 PSI = zeros(Ny,length(theta));
 
 switch M.type
  
  case {'arx','narx'}
   OO      = OPT; 
   OO.filt = 1;
   OO.fast = 1;
   for k=1:M.in(1).neta
    MM            = Mnew; 
    MM.in(1).type = 'linear';
    gg            = barx([y(:),z(:,k)],MM,OO);  
    PSI(:,k)      = -gg.phi*g.th;
   end
   index = length(g.th)+1;
   
   
  case {'fir','nfir'}
   OO      = OPT; 
   OO.filt = 1;
   OO.fast = 1;
   for k=1:M.in(1).neta,
    gg       = fir([y(:),z(:,k)],Mnew,OO);  
    PSI(:,k) = -gg.phi*g.th;
   end
   index = length(g.th)+1;
   
  case {'arma'}
   index = 1;    %  Keep track of where we are up to in filling up columns of PSI.

   % Derivatives w.r.t A
   for k=1:length(a)-1,
    num          = [zeros(1,k),1,zeros(1,length(a)-k-1)];
    psi          = ufilter(num,c,ygu,Mnew);
    PSI(:,index) = psi(:);
    index        = index+1;
   end
   
   % Derivatives w.r.t C
   for k=1:length(c)-1,
    num          = [zeros(1,k),1,zeros(1,length(c)-k-1)];
    psi          = -ufilter(num,c,pe,Mnew);
    PSI(:,index) = psi(:);
    index        = index+1;
   end
   
  case {'armax','narmax'}
   index = 1;    %  Keep track of where we are up to in filling up columns of PSI.
   
   % Derivatives w.r.t B
   for m = 1:nu  %  Loop over numerators for each input
    bit=(~strcmp(M.in(m).type,'linear') & M.nB(m)<1); % Just estimating a static nonlinearity?
    for k=0:M.nB(m)-bit
     num          = [zeros(1,k),1,zeros(1,M.nB(m)-k)];
     psi          = -ufilter(num,c,x(:,m),Mnew);
     PSI(:,index) = psi(:);
     index        = index+1;
    end
   end

   for k=1:M.nA
    num          = [zeros(1,k),1,zeros(1,M.nA-k)];
    psi          = ufilter(num,c,y,Mnew);
    PSI(:,index) = psi(:);
    index        = index+1;
   end

   % Derivatives w.r.t C
   for k=1:length(c)-1,
    num          = [zeros(1,k),1,zeros(1,length(c)-k-1)];
    psi          = -ufilter(num,c,pe,Mnew);
    PSI(:,index) = psi(:);
    index        = index+1;
            end
   
  otherwise %including nonlinear cases
   index = 1;    %  Keep track of where we are up to in filling up columns of PSI.
   
   % Derivatives w.r.t B
   if ~strcmp(lower(M.type),'arma'),
    for m = 1:nu  %  Loop over numerators for each input
     bit=(~strcmp(M.in(m).type,'linear') & M.nB(m)<1); % Just estimating a static nonlinearity?
     ac=conv(a(m,1:M.nA(m)+1),c);
     for k=0:M.nB(m)-bit
      num          = [zeros(1,k),1,zeros(1,M.nB(m)-k)];
      psi          = -w.*ufilter(conv(num,d),ac,x(:,m),Mnew);
      PSI(:,index) = psi(:);
      index        = index+1;
     end
    end

    for m = 1:nu  %  Loop over denominators for each input
     ac=conv(a(m,1:M.nA(m)+1),c);
     for k=1:M.nA(m)
      num          = [zeros(1,k),1,zeros(1,M.nA(m)-k)];
      psi          = w.*ufilter(conv(d,num),ac,pex(:,m),Mnew);
      PSI(:,index) = psi(:);
      index        = index+1;
     end
    end
   end

   % Derivatives w.r.t C
   for k=1:M.nC,
    num          = [zeros(1,k),1,zeros(1,length(c)-k-1)];
    psi          = -ufilter(num,c,pe,Mnew);
    PSI(:,index) = psi(:);
    index        = index+1;
   end

   % Derivatives w.r.t D
   for k=1:M.nD,
    num          = [zeros(1,k),1,zeros(1,length(d)-k-1)];
    psi          = ufilter(num,c,ygu,Mnew);
    PSI(:,index) = psi(:);
    index        = index+1;
   end
   
 end %END OF SWITCH FOR DIFFERENT MODEL TYPES

 % Now we handle gradient with respect to non-linearities
 zindex = 1; % Where we are up to in moving through columns of z
 for m=1:nu,  % Loop over all the inputs
  if any(strcmpi(M.type,{'ar','arma','arx','armax'})), nua = 1; else, nua = m; end
  for k=1:M.in(m).neta,
   psi          = ufilter(d,c,z(:,zindex),Mnew);
   psi          = ufilter(b(m,:),a(nua,:),psi,Mnew);
   PSI(:,index) = -psi(:);
   zindex       = zindex+1;
   index        = index+1;
  end
 end

 for k=1:M.out.neta,
  psi          = ufilter(d,c,zz(:,k),Mnew);
  PSI(:,index) = -psi(:);
  index        = index+1;
 end
 
 %Now prepare gradient and Jacobian terms
 phi  = PSI(OPT.n+1:Ny,:); 
 pe   = pe(OPT.n+1:Ny);
 grad = 2*phi'*pe/(Ny-OPT.n);

end;  % Test on whether div is set





