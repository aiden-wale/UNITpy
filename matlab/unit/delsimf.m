%  Delta operator version of dlsim.  That is, given a vector of inputs u() to
%  a plant expressed in delta operator transfer function form:
%
%  G(d) =     b_m d^m + ... b_1 d + b_0             d = q-1
%     	      -------------------------                -----
%      	      a_n d^n + ... a_1 d + a_0                delta
%
%  and with sampling period delta, work out the vector of outputs of the 
%  plant with zero initial conditions.  Usage is
%
%  y = delsimf(num,den,u,delta,y0)
%
%  where
%  
%  numd     = [b_m,...,b_0],  
%  dend     = [a_n,...,a_0].
%  delta    = sampling period in seconds
%  y0       = for n^th order system, y0 is a specification for the
%     	initial conditions y_0, y_1,...,y_{n-1}
%
%  This m file performs an identical function to delsim.m, but does it much
%  faster via the use of c code compiled to a mex function.
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

function y = delsimf(num,den,u,delta,y0)

if nargin<5 y0=[]; end;

MaxCount = length(u);
X = zeros(length(den) - 1,1);

%Put num, den into state-space form px=Ax+Bu (num and den are in descending powers of p)
[a,b,c,d] = tf2ss(num,den);

if (length(y0)>0)
 ad = eye(size(a))+delta*a;  bd = b*delta;
 T = c;  M = d*eye(size(a));
 for m=2:length(y0) T = [T;T(m-1,:)*ad]; end;
 h = T*bd;  h = [d;h(1:length(y0)-1)]; 
 M = toeplitz(h,[d,zeros(1,length(h)-1)]);
 uu = u(1:length(y0)); uu = uu(:);
 [U,S,V] = svd(T);
 jj = diag(S)>1e-12; jj = jj(:);
 S1 = inv(S(jj,jj)); V1T = V(:,jj); U1 = U(:,jj');
 x =  V1T*S1*U1'*(y0(:) - M*uu);
else
   x = 0*ones(size(b));
end;

% Pass to fast mex file engine to improve speed relative to Matlab
% interpreted loop
y = delsimeng(u,a,b,c,d,x,delta);







