% POLYNOM - function to pass a signal vector u through a nonlinearity X to
% produce an output x = X(u) where X is defined as a polynomial.  That
% is, a vector alpha = [a1,a2,...,an] is specified and then
%
% X = a1*u + a2*u^2 + a3*u^3 + ... + an*u^n
%
%
% Usage is: 
%
% [x,phi] = polynom(u,alpha)
%
% Where
%
% u     = input vector
% alpha = parameters vector specifying X.
% x     = vector representing output of nonlinearity
% phi   = matrix of `regressors' in the sense that x = phi*alpha(:)
% w     = derivative dX(u)/du
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

function [x,phi,w] = polynom(u,alpha)

u = u(:); 
wun = ones(1,length(alpha)); 
index = 0:1:length(alpha)-1;
phi = (u*wun).^kron(index,ones(length(u),1));
x = phi*alpha(:);
beta = alpha(:); 
beta = alpha(2:length(alpha)); 
beta = beta(:); kk = 2:1:length(alpha); 
beta = beta.*kk(:);
w = alpha(1)*ones(size(u))+phi(:,1:length(alpha)-1)*beta;
