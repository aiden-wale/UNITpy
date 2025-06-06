%  SSTOTF: Function to add the polynomial form descriptions to a model
%  structure that are equal to the polynomial form descriptions given in
%  the model stucture 
%
%  Usage is 
%
%  g = sstotf(G);
% 
%  where
%
%   G         = Initial model structure specification, which should
%               contain elements G.ss.A, G.ss.B, G.ss.C, G.ss.D and
%               G.ss.K that specify the innovations form model
%
%               x_{t+1} = Ax_t + Bu_t + Ke_t
%               y_t = Cx_t + Du_t + e_t
%
%   g         = Given model structure G, with elements G.A, G.B,G.C, G.D
%               added/augmented/changed as need be so that thy represent
%               the polynomial form description
%
%               y_t = B/A u_t + C/D e_t
%
%               which is steady state input-output equivalent to the
%               state space system specified in the input data structure.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                     		      Australia.

% Copyright (C) Brett Ninness.

function g=sstotf(G);

% Copy all params in input structure to output one
g = G;

% Increase robustness to poorly specified input
if ~isfield(G.ss,'K') G.ss.K=[]; end;
if ~isfield(G.ss,'D') G.ss.D=[]; end;

% Now overwrite the transfer function bits according to ss bits
if isfield(G,'ss')
 [nx,nu]=size(G.ss.B); [ny,nx]=size(G.ss.C); 
 g.A = zeros(nu,nx+1,ny); g.B=g.A; 
 if isempty(G.ss.D); G.ss.D=zeros(ny,nu); end;
 for k=1:nu   % Dynamics model
  for m=1:ny
   [g.B(k,:,m),g.A(k,:,m)] = ss2tf(G.ss.A,G.ss.B,G.ss.C(m,:),G.ss.D(m,:),k);
  end;
 end;
 if ~isempty(G.ss.K)
  g.C = zeros(ny,nx+1,ny); g.D=g.C;  
  for k=1:ny  % Noise model
   for m=1:ny
    [g.C(k,:,m),g.D(k,:,m)] = ss2tf(G.ss.A,G.ss.K,G.ss.C(m,:),eye(1,ny),k);
   end;
  end;
 else
  g.C = zeros(ny,1,ny); g.D=g.C;  
  for k=1:ny 
   for m=1:ny g.C(k,1,m)=1; g.D(k,1,m)=1; end; end;
 end;
 % Finally, fill in bits specifying orders of numerators and denominators
 g.nA=nx*ones(nu,1); 
 if ~strcmpi(G.par,'grey')  % With grey box par, nB ~= nA is possible
  g.nB=g.nA;
 end; 
else
 error('Need to supply ss model to sstotf')
end;