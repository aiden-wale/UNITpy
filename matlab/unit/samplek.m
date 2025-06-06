%  SAMPLEK: This routine takes a continuous time stochastic state 
%  space model description
%
%  dx(t) = Ax(t)dt + Bu(t)dt + Ke(t) 
%  dy(t) = Cx(t)dt + e(t)
%
%  Cov(e(t)) = I, x(0) = X1.
%
%  and, under the assumption of zero order hold sampling on u(t), and 
%  integrated sampling of y(t) over the integration period T equal to the
%  sampling period T, delivers the discrete time stochastic state 
%  space description
%
%  x_{k+1} = Fx_k + Gu_k + Lw_k
%      y_k = Hx_k + Du_k + v_k
%
%  Cov(w_k) = Q; Cov(v_k) = R; Cov(w_kv_k^T) = S;
%
%  using the methods presented in the paper [1]:
%
%  [1] "Issues in sampling and estimating continuous time models with
%       stochastic disturbances"  Lennart Ljung and Adrian Wills.
%       Automatica, Vol. 46, No. 5, pp. 925 - 931, 2010.
%
%  Usage is:
%
%  M = samplek(m);
%
%  Where
%
%  m      = model structure specifying continuous time state space
%           system via m.ss.A, m.ss.B, m.ss.C, m.ss.K, m.ss.X1, m.T
%
%  M      = Discrete time model valid for t_k < t < t_k+m.T where
%  M.ss.A = F
%  M.ss.B = G
%  M.ss.C = H
%  M.ss.D = D
%  M.ss.K = L
%  M.ss.Q = Q
%  M.ss.R = R
%  M.ss.S = S
%
%  Written by Brett Ninness,  School of EE & CS
%             Adrian Wills    University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness

function M = samplek(mm)

% Extract continuous time system matrices and sampling period.

T   = mm.T;
D   = T;
a   = mm.ss.A;
b   = mm.ss.B;
c   = mm.ss.C;
if isempty(mm.ss.K),
 r1  = zeros(size(a));
 r12 = zeros(size(a,1),size(b,2));
else
 r1  = mm.ss.K*mm.ss.K';
 r12 = mm.ss.K;
end
r2  = eye(size(mm.ss.C,1));
if isempty(mm.ss.X1),
 x1 = zeros(size(a,1),1);
else
 x1  = mm.ss.X1;
end
 

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

try
 % Compute Kalman gain from Riccati equation
 [X,L,KK]=dare(f',h',(q+q')/2,(r+r')/2,s,i);
 Knew = KK';
catch err
 Knew = zeros(size(c'));
end

M       = mm;
M.ss.A  = f;
M.ss.B  = g;
M.ss.C  = h;
M.ss.D  = d;
M.ss.K  = Knew;
M.ss.Q  = q;
M.ss.R  = r;
M.ss.S  = s;
M.ss.X1 = x1;
