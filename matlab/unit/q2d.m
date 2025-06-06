%  Function to convert forward shift operator discrete time 
%  model to a delta operator equivalent.  
% 
%  Usage is:   [n,d] = q2d(num,den,T);
%
%  where
% 
%  num,den = shift operator system in transfer function form
%        T = sampling period in seconds
%      n,d = delta operator equivalent sytsem in tf form
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

function [n,d] = q2d(num,den,T);

[a,b,c,d] = tf2ss(num(:)',den(:)');
a = (a - eye(size(a)))/T;
b = b/T;
[n,d] = ss2tf(a,b,c,d,1);

