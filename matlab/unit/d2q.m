%  Function to convert delta operator discrete time 
%  model to a shift operator equivalent.  Philosophically
%  there is not reason to do this, and it is inherently a bad thing to do
%  on numerical grounds at fast sampling rates.  Nevertheless, sometimes
%  it's useful in comparison excercises.
% 
%  Usage is:   [n,d] = d2q(num,den,T);
%
%  where
% 
%  num,den = delta operator system in transfer function form
%        T = sampling period in seconds
%      n,d = shift operator equivalent sytsem in tf form
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness 

function [n,d] = d2q(num,den,T);

[a,b,c,d] = tf2ss(num(:)',den(:)');
a = a*T + eye(size(a));
b = b*T;
[n,d] = ss2tf(a,b,c,d,1);