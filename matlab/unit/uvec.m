% This is a function that returns a unit column vector of length n with
% entry 1<=k<=n equal to 1.0 and all other entries equal to 0.
%
% Usage:
%
%   v = uvec(k,n);
%
% Written by Adrian Wills: School of EE & CS
%                          University of Newcastle
%         		           Australia.

% Copyright (C) Adrian Wills

function v = uvec(k,n)

v = zeros(n,1); v(k) = 1;