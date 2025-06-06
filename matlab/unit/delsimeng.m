% Function to perform filtering according to delta operator state space 
% description of filter.  
%
%    Written by Brett Ninness, School of Electrical Engineering
%                              and Computer Science
%                              University of Newcastle
%                              Australia.

% Copyright (C) Brett Ninness.


function y=delsimeng(u,a,b,c,d,x,delta)

%Simulate a state-space system given in delta operator form...

%Number of data points, inputs, outputs and states

N  = max(size(u));

if isempty(a), a = zeros(length(x)); end
if isempty(b), b = zeros(length(x),min(size(u))); end
if isempty(c), c = zeros(1,length(x)); end
if isempty(d), d = zeros(1,min(size(u))); end

nu = min(size(u));
ny = size(c,1);
nx = size(a,1);

a = delta*a; 
b = delta*b;

u = u.';
y = zeros(ny,N);
for i = 1:N;
 y(:,i) = y(:,i) + d*u(:,i);
 if nx>0,
  y(:,i) = y(:,i) + c*x;
  x_increment = a*x + b*u(:,i);
  x = x + x_increment;
 end
end

y = y.';