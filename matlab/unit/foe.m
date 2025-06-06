%  Computes an estimate using an output-error model
%
%  F(w_s) = B(p)/A(p)cos(2*pi*w_s) + V(w_s)
%
%  from frequency domain data observations F(w_s) that are corrupted by
%  additive noise V(w_s).  The operator p can be the Z tranform variable
%  with z=e^(j*w_s*T) (with T being the sampling period in seconds), the
%  Euler differencing (delta) operator d = (q-1)/T with d =
%  (e^(j*w_s*T)-1)/T or the Laplace Transform variable s with s=j*w.  A
%  quadratic (least squares) loss criterion is used.
%
%  Usage is 
%
%  G = foe(Z,M,OPT);
% 
%  Z          = observed frequency response data [F(:),w(:)] where
%               plot(w,abs(F)) should plot the measured  frequency
%               response.   Units for w are real *not* normalised freq.
%  M          = Data structure which defines the model structure which
%               is to be estimated from the data as follows:
%   M.A       = Number of poles to be estimated in denominator - which is
%               then set as equal to # of zeros to be estimated in numerator.
%   M.op      = set to 'q' for shift, 'd' for delta, 's' for Laplace  
%               Default = 's'.
%   M.T       = sampling period in s. (Ignored for q case) Default = 1;
%
%  OPT        = Data structure which defines options for the estimation
%               algorithm as follows:
%   OPT.basis = only applicable for 's' operator models, and selects either
%               Chebychev ('cheby') or Laguerre ('ortho') orthonormal bases,
%               or normal non-orthonormal polynomial ('polyb') basis.  
%               Default is 'ortho'.
%   OPT.W     = Vector of same dimension as w that specifies a
%               frequency weighting for the least squares fit.  That is,
%               plot(w,W) should give a graphical interpretation of the
%               weighting.  The default is a flat (unprejudiced)
%               weighting.
%   OPT.mdec  = Minimum relative decrease of cost before search is
%               terminated.  Default is 1e-8;spl
%
%  G          = Data structure which specifies the estimated model as
%               follows:
%   G.B/G.A   = estimated transfer function for model of dynamics.
%   G.G       = Frequency response of estimated model for dynamics.
%   G.th      = Estimated Parameter vector from which G.B, G.A are formed. 
%
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness.

function G = foe(Z,M,OPT)

% Extract out relevant vectors from input data
[F,w,ny,nu,Ny] = Z2data(Z); F=squeeze(F); F=F(:); wmax = max(w);

% Check what options not specified explicitly by user and set to defaults

if ~exist('OPT') 
  OPT.dsp=0; OPT.miter=20; OPT.tol=1e-5; 
  OPT.lmax=20; OPT.basis = 'ortho'; OPT.mdec = 1e-8;
  OPT.W = ones(size(F));  
else
  if ~isfield(OPT,'dsp')   OPT.dsp=0;           end;
  if ~isfield(OPT,'miter') OPT.miter=20;        end;  
  if ~isfield(OPT,'tol')   OPT.tol=1e-5;        end;    
  if ~isfield(OPT,'lmax')  OPT.lmax=20;         end;      
  if ~isfield(OPT,'mdec')  OPT.mdec=1e-8;       end;        
  if ~isfield(OPT,'basis') OPT.basis='ortho';   end;
  if ~isfield(OPT,'W')     OPT.W=ones(size(F)); end;   
  if (length(OPT.basis) ~=5) 
    error('Not a recognised basis from: ortho,cheby,polyb'); end;      
end;

% Check which parts of model structure were unspecified and set to defaults.

if ~exist('M') error('Need to specify initial model structure M!'); 
  else
  if ~isfield(M,'op')    M.op='q';   end;
  if ~isfield(M,'T')     M.T=1;      end;  
  if ~isfield(M,'B')     M.B=M.A;    end;  
  if ~isfield(M,'delay') M.delay=0;  end;    
  if ~isfield(M,'w')     M.w= logspace(log10(pi/M.T/1000),log10(pi/M.T)); end;      
  M.A = M.A(:); M.B = M.B(:);
end;

%  Establish frequency domain variable appropriate to time domain operator
if (M.op=='q') ww = exp(j*M.w*M.T); 
elseif (M.op=='d') ww = (exp(j*M.w*M.T)-ones(size(M.w)))/M.T; 
else ww = j*M.w; end;  

%  Is frequency normalisation necessary?
normw=0;  if [M.op == 's', OPT.basis ~= 'ortho'] normw = 1; end;

% Check to see of only integer orders where specified as initial guesses
% for dynamics: if so get initial estimate by fitting ARX model structure.

if [length(M.A)<2  floor(M.A(:)')==M.A(:)']  
 % Get initial ARX estimate of specified order;      
 g = farx(Z,M,OPT); th0 = g.th; X = g.X; 
 M.A = g.A; M.B = g.B; n = length(g.B);
else  % Otherwise, re-express initial guess wrt chosen basis
 ff = polyval(M.B,ww)./polyval(M.A,ww);  % Response of initial guess
 % Use farx to translate this initial guess to requested basis.
 g = farx([ff(:),M.w(:)],M,OPT); th0 = g.th; X = g.X;  
 M.A = g.A; M.B = g.B; n = length(g.B); 
end;

%  Now use iterative Gauss-Newton search to find minimum of quadratic cost.

th = argmin(Z,'VNf',th0,OPT,M);

%  Extract A and B polynomials from theta parameter vector

G.B = th(1:n)'; G.A = th(n+1:length(th))'; 

%  Revert basis to normal polynomial one.

G.B = G.B*X; G.A = [1,G.A]*X;

%  Undo frequency normalisation if it was applied

if (normw)
  G.A = G.A.*(wmax.^(1:n));   %  Undo frequency Normalisation if 
  G.B = G.B.*(wmax.^(1:n));   %  it was applied
end;  

% Pack results into output data structure.

G.B = G.B/G.A(1); G.A = G.A/G.A(1); G.delay=M.delay; 
G.T = M.T; G.w = M.w; G.op = M.op; G.th = th; G.type='foe';
G.C=[]; G.D=[];

% Add legend for prospective plotting
G.disp.legend=['Estimated ',G.type,' model'];

G.alg='gn'; % Record that Gauss-Newton search was employed

















