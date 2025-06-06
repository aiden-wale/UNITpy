%    Delta operator version of filter.  That is, given a vector of inputs u() to
%    a plant expressed in delta operator transfer function form:
% 
%    G(d) =     b_0+b_1 d^{-1}+...+b_m d^{-m}             d = q-1
%               -----------------------------                -----
%               a_0+a_1 d^{-1}+...+a_n d^{-n}                delta
% 
%    and with sampling period delta, work out the vector of outputs of the 
%    plant with zero initial conditions.  Usage is
% 
%    y = delfilter(num,den,u,delta,y0)
% 
%    where
%    
%    numd     = [b_0,...,b_m],  
%    dend     = [a_0,...,a_n].
%    delta    = sampling period in seconds
%    y0       = for n^th order system, y0 is a specification for the
%               initial conditions y_0, y_1,...,y_{n-1}
%
%    Written by Brett Ninness, School of Electrical Engineering
%                              and Computer Science
%                              University of Newcastle
%                              Australia.

% Copyright (C) Brett Ninness.

function y = delfilter(num,den,u,delta,y0)

% Setup numerators and denominators in d form that is equivalent to d^{-1} form specified
if nargin<5 y0=[]; end;
dabk = length(num)-length(den); % delfilt in d^{-1}, delsimf in d => zero padding

% NOTE: for large relative degrees then poly(-ones(dabk,1)) or
% poly(-ones(-dabk,1)) will have very large coefficients - a problem with
% delta operator.
b = conv(num(:)',poly(-ones(-dabk,1))); 
a = conv(den(:)',poly(-ones(dabk,1))); 

% OK, we've transformed from d^{-1} form to equivalent d form - use code
% for that form
y = delsimf(b,a,u,delta,y0);










