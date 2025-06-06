%   Calculates the chi squared distribution function, as opposed to 
%   the chi squred probability density function which is in chi.m
%
%   Usage:	y = chidist(x,d);
%
%   where
% 
%   x:		X axis argument
%   d:		Degrees of freedom of the chi squared distribution function
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness

function y = chidist(x,d)

if x > 0.0
 y = quad8('chi',0,x,[],[],d);
else
 y = 0.0;
end;
