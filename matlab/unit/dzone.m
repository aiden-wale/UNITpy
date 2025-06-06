%  Function to deadzone a signal u between a 
%  lower and upper limit
%
%  Usage is function z = dzone(u,lower,upper)
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

function z = dzone(u,lower,upper)

[m,n]=size(u); u=u(:); z=zeros(size(u));
z(u<lower) = u(u<lower)-lower;
z(u>upper) = u(u>upper)-upper;
if n>m, z=z'; end




