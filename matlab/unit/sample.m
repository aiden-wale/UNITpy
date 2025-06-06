%   SAMPLEK: This routine takes a continuous time state space description
%
%   dx(t) = Ax(t)dt + Bu(t)dt + w(t) 
%   dy(t) = Cx(t)dt + v(t)
%
%   and, under the assumption of zero order hold sampling on u(t), and 
%   integrated sampling of y(t) over the integration period m.T, delivers
%   the discrete time state space model
%
%   x_{k+1} = Fx_k + Gu_k + v_k
%       y_k = Hx_k + Du_k + z_k
%
%   using the methods presented in the paper [1]:
%
%   [1] "Issues in sampling and estimating continuous time models with
%        stochastic disturbances"  Lennart Ljung and Adrian Wills.
%        Automatica, Vol. 46, No. 5, pp. 925 - 931, 2010.
%
%   Usage is:
%
%   M = samplek(m);
%
%   Where
%
%   m      = model structure specifying continuous time state space
%            system via m.ss.A, m.ss.B, m.ss.C, m.ss.K, m.ss.X1, m.T
%
%   M      = Discrete time model valid for t_k < t < t_k+m.T where
%   M.ss.A = F
%   M.ss.B = G
%   M.ss.C = H
%   M.ss.D = D
%
% Written by Brett Ninness,  School of EE & CS
%            Adrian Wills    University of Newcastle
%             		         Australia.

% Copyright (C) Brett Ninness

function [f,g,h,d,q,s,r] = sample(M,theta,T,D)

if nargin<4,
	D=T;
end

% Make call to pf to get initial matrices
mm  = theta2m(theta,M,1);
a   = mm.ss.A;
b   = mm.ss.B;
c   = mm.ss.C;
r1  = mm.ss.Q;
r12 = mm.ss.S;
r2  = mm.ss.R;

% Get number of states and construct big matrix in a similar manner to [1]
n  = size(a,1);
i  = eye(n);
z  = zeros(n);
C  = [-a   i   z   z   z;
    z  -a  r1   z   z;
    z   z  a'   i   z;
    z   z   z   z   i;
    z   z   z   z   z];
   
% Take matrix exponential and extract key sub-matrices   
eC = expm(C*D);
f3 = eC(2*n+1:3*n,2*n+1:3*n);
g2 = eC(n+1:2*n,2*n+1:3*n);
g3 = eC(2*n+1:3*n,3*n+1:3*n+n);
h2 = eC(n+1:2*n,3*n+1:3*n+n);
h3 = eC(2*n+1:3*n,4*n+1:5*n);
k1 = eC(1:n,3*n+1:4*n);

% Form part of sampled system that depends on integration time d
s  = (f3'*h2*c'+g3'*r12)/D;
h  = (c*g3')/D;
d  = (c*h3'*b)/D;
r  = (c*(f3'*k1 + (f3'*k1)')*c' + c*h3'*r12 + r12'*h3*c')/D/D + r2/D;

% Check if T==D (roughly)
if abs(T-D)>100*eps,
 % If not then compute new version of the above but using T instead of d
 eC = expm(C(1:4*n,1:4*n)*T);
 f3 = eC(2*n+1:3*n,2*n+1:3*n);
 g2 = eC(n+1:2*n,2*n+1:3*n);
 g3 = eC(2*n+1:3*n,3*n+1:3*n+n);
 h2 = eC(n+1:2*n,3*n+1:3*n+n);
 f  = f3';
 g  = g3'*b;
 q  = f3'*g2;
else
 % T==d so we have all the matrices we need to compute state transition
 f  = f3';
 g  = g3'*b;
 q  = f3'*g2;
end
