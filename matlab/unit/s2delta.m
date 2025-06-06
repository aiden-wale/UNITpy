% 	Converts continuous time transfer function to delta domain
%	discrete time equivalent assuming zero order hold on input.
%
% [numd,dend] = s2delta(num,den,deltat)
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

function [numd,dend] = s2delta(num,den,deltat)
[a,b,c,d] = tf2ss(num,den);

if deltat ~= 0                      % If delta=0 then no change occurs. 
 [m,n] = size(a);
 [m,nb] = size(b);
 I = eye(n,n); O = 0*I;
 omega = [I O]*expm(([a I; O O])*deltat)*[O;I]/deltat;
 f = omega*a;
 g = omega*b;
end

h = c; k = d; [numd,dend] = ss2tf(f,g,h,k,1);

