%  SSTOTF: Function to add the polynomial form descriptions to a model
%  structure that are equivalent to the existing state space form form
%  descriptions.
%
%  Usage is
%
%  g = sstotf(G);
%
%  where
%
%   G = Initial model structure specification, which should
%       contain elements G.ss.A, G.ss.B, G.ss.C, G.ss.D and
%       G.ss.K that specify the innovations form model
%
%       x_{t+1} = Ax_t + Bu_t + Ke_t
%       y_t = Cx_t + Du_t + e_t
%
%   g = Given model structure G, with elements G.A, G.B,G.C, G.D
%       added/augmented/changed as need be so that they represent
%       the polynomial form description
%
%       y_t = B/A u_t + C/D e_t
%
%       which is steady state input-output equivalent to the
%       state space system specified in the input data structure.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                         Australia.

% Copyright (C) Brett Ninness.

function g=tftoss(G);

% Copy all params in input structure to output one
g = G;

% Now overwrite the ss bits according to the transfer function bits
% Dynamics model first
if [isfield(G,'A'), isfield(G,'B')]
 ny=size(G.A,3);  % Find out how many outputs we have
 nu=size(G.B,1);  % Find out how many inputs we have
 % Now determine maximum state dim - complicated by fact of order(B)>order(A)?
 ndim=0;
 for m=1:ny
  for k=1:nu 
   ndim = ndim+max(size(G.B,2)-1,size(G.A,2)-1); 
  end
 end;
 A = zeros(ndim,ndim);  % Start out with maximum dim matrices
 B = zeros(ndim,nu);
 C = zeros(ny,ndim);
 D = zeros(ny,nu); xidx = 1;
 for m=1:ny  % Fill in blocks for each i/o pair
  for k=1:nu
   if any(strcmpi(G.type,{'ar','arma','arx','armax'})), nua = 1; else, nua = k; end
   % First compute any padding necessary of more lags in num than den
   zpad = zeros(1,length(G.B(k,:,m))-length(G.A(nua,:,m)));
   [a,b,c,d]=tf2ss(G.B(k,:,m),[G.A(nua,:,m),zpad]);
   [nax,nax]=size(a);  % Find state dimension of this i/o pair;
   A(xidx:xidx+nax-1,xidx:xidx+nax-1) = a;  % Put it into augmented system
   B(xidx:xidx+nax-1,k) = b;
   C(m,xidx:xidx+nax-1) = c;
   D(m,k) = d;
   xidx = xidx+nax;  % Update record of where next block should go
  end;
 end;
 % Now need to cut things down to minimal state dimension
 A = A(1:xidx-1,1:xidx-1);  B = B(1:xidx-1,:); C = C(:,1:xidx-1);
 %[g.ss.A,g.ss.B,g.ss.C,g.ss.D] = minreal(A,B,C,D);
 [g.ss.A,g.ss.B,g.ss.C,g.ss.D] = modred(A,B,C,D);
else
 g.ss.A=[]; g.ss.B=[]; g.ss.C=[]; g.ss.D=[];
end;

if any(strcmpi(G.type,{'ar','arma','arx','armax'})), G.D = G.A; end

% Then the noise model
if [isfield(G,'C'), isfield(G,'D')]
 if [~isempty(G.C) ~isempty(G.D)]
  nx = size(G.A,2);         % Maximum number of poles on any input model
  ndim = (nx-1)*ny*ny;      % Maximum possible noise state dimension
  A = zeros(ndim,ndim);     % Start out with maximum dim matrices
  K = zeros(ndim,ny);
  C = zeros(ny,ndim);
  D = zeros(ny,ny); xidx = 1;
  for m=1:ny  % Fill in blocks for each i/o pair
   for k=1:ny
    [a,b,c,d]=tf2ss(G.C(k,:,m),G.D(k,:,m));
    [nax,nax]=size(a);  % Find state dimension of this i/o pair;
    A(xidx:xidx+nax-1,xidx:xidx+nax-1) = a;  % Put it into augmented system
    K(xidx:xidx+nax-1,k) = b;
    C(m,xidx:xidx+nax-1) = c;
    D(k,k) = 1;         % Chance of a bug creeping in here if a G.C or G.D not monic!
    xidx = xidx+nax;    % Update record of where next block should go
   end;
  end;
  % Now need to cut things down to a minimal state dimension
  A = A(1:xidx-1,1:xidx-1);  K = K(1:xidx-1,:); C = C(:,1:xidx-1);

  % Augment this to any dynamics model
  z1 = zeros(size(g.ss.A,1),size(A,2));  z2 = zeros(size(A,1),size(g.ss.A,2));
  z3 = zeros(size(g.ss.B,1),size(K,2));  z4 = zeros(size(K,1),size(g.ss.B,2));

  % Grab out bits just to do with i/o dynamics
  aa = [g.ss.A,z1;z2,A]; bb = [g.ss.B,z3;z4,K]; cc = [g.ss.C,C]; dd = [g.ss.D,D];
  %[g.ss.A,b,g.ss.C,d]=minreal(aa,bb,cc,dd); 
  [g.ss.A,b,g.ss.C,d]=modred(aa,bb,cc,dd); 
  g.ss.B = b(:,1:nu); g.ss.D = d(1:ny,1:nu);

  % Grab Kalman gain
  g.ss.K=b(:,nu+1:end);
 else % No noise model => zero noise spec in ss domain
  g.ss.K=[];  g.ss.R = zeros(size(ny,ny),1);
 end;
else
 g.ss.K=[];
end;

