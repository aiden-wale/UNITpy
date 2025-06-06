%   FSID computes a state space model [A,B,C,D] that fits an observed (possibly
%   multi-output) frequency response F.  This routine is an implementation
%   of the algorithm developed by McKelvey, Akcay and Ljung, IEEE
%   Transactions on Automatic Control, V41(7), pp960-979, 1996.  See also
%   paper by same authors, Automatica V32(6), pp885-902, 1996.
%
%   The estimated model [A,B,C,D] can be found in either discrete time
%   shift operator form, or continuous time form.
%
%   Usage is:  G = fsid(Z,M,OPT);
%
%   where
%
%   Z         = observed frequency response data [F(:),w(:)] where
%               plot(w(:,1),abs(F)) should plot the measured  frequency
%               response.   Units for w are real *not* normalised freq,
%               measured in radians per second.
%   M         = Definiton of the model structure which
%               is to be estimated from the data as follows:
%    M.A      = Number of poles to be estimated in denominator - which is
%               then set as equal to # of zeros to be estimated in numerator.
%    M.op     = set to 'q' for shift or 's' for Laplace.  Default is M.op='q';
%    M.T      = sampling period in s. (Ignored for q case) Default = pi/max(w);
%   OPT       = Data structure which defines options for the estimation
%               algorithm as follows:
%    OPT.lag  = `Embedding' dimension for state-space - dimension of
%               space in which to search for M.A'th dimensional subspace
%               that the state lives in.  Default is OPT.lag = 4*M.A;
%    OPT.R    = (Optional) vector of noise variances at each frequency
%               measurement indexed in w.  Default is OPT.R set to all 1's so
%               that all measurements are weighted equally.
%  G          = Data structure which specifies the estimated model as
%               follows:
%   G.B,G.A   = State Space representation of estimated model
%   G.C,G.D
%   G.G       = Frequency response of estimated model - a column of G.G
%               is generated corresponding to each output column
%               specified in Z.
%   G.sing    = Singular values that arise in calculating rank of systems
%               observability matrix - examining this can give an
%               indication of the underlying system order.
%
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%        		              Australia.

% Copyright (C) Brett Ninness.

function [g] = fsid(z,mm,OPT)

% Make sure data is OK
z=startZ(z);

% Now extract things from data structure
% Get sizes of stuff: p=# outputs, m=# inputs, M=# frequency points
[G_in,w,p,m,M] = Z2data(z);

% Check which parts of model structure were unspecified and set to defaults.
if ~exist('mm'),
 error('Need to specify initial model structure M!');
elseif ~isfield(mm,'nx'),
 error('Need to at least specify M.nx');
elseif isempty(mm.nx),
 error('M.nx cannot be empty.');
else
 if ~isfield(mm,'op'),    mm.op='s';                                          end;
 if ~isfield(mm,'T'),     mm.T = pi/max(w);                                   end;
 if ~isfield(mm,'delay'), mm.delay=zeros(m,1);                                end;
 if ~isfield(mm,'w'),     mm.w= logspace(log10(pi/mm.T/1000),log10(pi/mm.T)); end;
end;

n = mm.nx;

if ~exist('OPT')  % Default horizon is four times selected order
 OPT.lag = min(4*n,floor(M/2));
 OPT.R   = ones(1,length(w));    %  Default noise weighting is none.
else
 if ~isfield(OPT,'lag'), 
  OPT.lag = min(4*n,floor(M/2)); 
 end
 if ~isfield(OPT,'R'),
  OPT.R = ones(1,length(w));
 else
  if (length(OPT.R) ~= length(w))
   error('Lengths of M.w and OPT.R must match');
  end
 end
end

if (OPT.lag<n)
 error('Must have model order less than OPT.lag')
end;

% Make variables local
i = OPT.lag;

% Perform the pre-warp if continuous
if strcmp(mm.op,'s')
 T  = 2*pi/max(w); 
 Tq = 1; 
 w  = 2*atan(w*T/2);
else
 Tq = mm.T;
end

% Generate exp(j*w*Tq) and G and Wm
ew = exp(j*Tq*w); 
ew = ew(:).'; 
G  = zeros(i*p,M*m);
for r=1:i, 
 for k=1:M, 
  G((r-1)*p+1:r*p,(k-1)*m+1:k*m) = exp(j*Tq*w(k)*(r-1))*G_in(:,:,k); 
 end
end
W  = power(kronaw(ones(i,1),ew),kronaw(ones(1,M),[0:1:i-1]')); 
Wm = kronaw(W,eye(m));
Wp = kronaw(W,eye(p));
if any(OPT.R~=1)  % Check to see if non-trivial weighting specified
 WRp=zeros(size(Wp));
 for k=1:length(OPT.R), 
  WRp(:,k)=sqrt(OPT.R(k))*Wp(:,k); 
 end
 K = triu(qr(WRp')); 
 K = real(K(1:i,1:i))';
else
 K = eye(p*i,p*i);
end;

%  Perform projection by QR factorization
R   = triu(qr([real(Wm) imag(Wm);real(G) imag(G)].')).';
R22 = full(R(m*i+1:m*i+p*i,m*i+1:end))';

% Form SVD
[U,S,V] = svd(R22/K,0); 
Us = V(:,1:n);

% Construct A and C
A = pinv(Us(1:(i-1)*p,:))*Us(p+1:end,:); 
C = Us(1:p,:);

% Construct B and D
CX  = frmimo(A,eye(n),C,exp(j*w*Tq)); 
PHI = zeros(M*p,n+p); 
RHS = zeros(M*p,m);
for k=1:M,
 PHI((k-1)*p+1:k*p,:) = [CX(:,:,k) , eye(p)];
 RHS((k-1)*p+1:k*p,:) = G_in(:,:,k);
end
BD = pinv([real(PHI);imag(PHI)])*[real(RHS);imag(RHS)];
B  = BD(1:n,:); 
D  = BD(n+1:end,:);

% If continuous then reconstruct the analog system (a,b,c,d)
g=mm;
if strcmp(mm.op,'s'),
 aa = inv(eye(size(A))+A);
 a  = (2/T) * aa * (A-eye(size(A)));
 b  = (2/sqrt(T))*aa*B;
 c  = (2/sqrt(T))*C*aa;
 d  = D - C*aa*B;
 g.ss.sys=ss(a,b,c,d);
else
 a=A; b=B; c=C; d=D;
 g.ss.sys=ss(a,b,c,d,mm.T);
end


%Now make sure g has relevant structure entries
g.ss.A  = a;
g.ss.B  = b; 
g.ss.C  = c; 
g.ss.D  = d;
g.ss.K  = [];
g.ss.F  = [];
g.ss.G  = [];
g.ss.X1 = [];
g.T     = mm.T; 
g.w     = mm.w; 
g.op    = mm.op; 
g.sing  = diag(S); 
g.delay = zeros(m,1);
g.type  = 'ss';
g.par   = 'full';
g       = sstotf(g);
g.var   = 1;


% Add legend for prospective plotting
g.disp.legend=['Estimated ',upper(g.type),' model via Subspace Alg.'];

g.alg='sid'; % Record that block solution was used