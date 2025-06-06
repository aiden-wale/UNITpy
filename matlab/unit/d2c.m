%  D2C: Function to convert discrete time model structure to continuous
%  time one.
%
%  Usage is
%
%  G = d2c(M);
%
%  where
%
%  M = Discrete time model structure
%
%  G = Continuous time structure assuming sampling
%      rate of M.T seconds.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                     		  Australia.

% Copyright (C) Brett Ninness.

function G=d2c(M)

G = M;   % Copy all input elements over

% Convert d->c via inverse Tustin transform
A = M.ss.A; B = M.ss.B; C = M.ss.C; D = M.ss.D;
aa      = inv(eye(size(A))+A);
G.ss.A  = (2/M.T) * aa * (A-eye(size(A)));
G.ss.B  = (2/sqrt(M.T))*aa*B;
G.ss.C  = (2/sqrt(M.T))*C*aa;
G.ss.D  = D - C*aa*B;

% Output the d2c result in TF form as well
[G.B,G.A] = ss2tf(G.ss.A,G.ss.B,G.ss.C,G.ss.D,1);

% New operator is continuous time one
G.op = 's';
