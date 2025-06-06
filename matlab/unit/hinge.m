% HINGE - function to pass a signal vector u through a nonlinearity X to
% produce an output x = X(u) where X is defined by a set of `hinging
% hyperplanes'.  That is, a vector alpha made of up groups of pairs of
% parameters a11, a12, a21, a22,....,am1,am2 specifies X as
%
% X = a11+a12*u + sum_{k=2}^m X_k(u,ak1,ak2)
%
% where  X_k(u,ak1,ak2) = ak1 + ak2*u ; u >  -ak1/ak2
%                       = 0           ; u <= -ak1/ak2
%
% Usage is:
%
% [x,w] = hinge(u,alpha)
%
% Where
%
% u     = input vector
% alpha = parameters vector specifying X.
% x     = vector representing output of nonlinearity
% w     = derivative dX(u)/du
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

function [x,w] = hinge(u,alpha)

m = length(alpha);
if rem(m,2) error('alpha must contain pairs of parameters'); end;

wun = ones(size(u)); 
x   = alpha(1)*wun + alpha(2)*u; 
w   = alpha(2)*ones(size(u));
for k = 1:m/2-1
 a1 = alpha(2*k+1); 
 a2 = alpha(2*k+2);
 if (abs(a2)>eps), 
 	breakpoint = -a1/a2; 
 else
 	breakpoint = 1e20; 
 end
 index       = logical(u > breakpoint*wun);
 xnew        = zeros(size(u));
 wnew        = xnew;
 xnew(index) = a1*wun(index) + a2*u(index);
 wnew(index) = a2*wun(index);
 x           = x+xnew; 
 w           = w+wnew;
end;






