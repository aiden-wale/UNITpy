function z = sat(u,lower,upper,m)
% Function to saturate a signal u in the 
% interval [lower,upper] and with gain m.
%
%  Usage is function z = sat(u,lower,upper,m)
%
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

if nargin<3 m=1; end;
if nargin<2 upper=-lower; end;
wun = ones(size(u));
z = max(min(m*u,upper*wun),lower*wun);

