%  Function to compute frequency response of LTI system in state space form
%
%  Usage G = mimofr(A,B,C,s)
%
%  In which case
%
%   G = C*[sI-A)^-1]*B
%
%  for each of the M entries in s.  If B and C imply ny outputs and nu
%  inputs (ie, a multivariable (MIMO) model) then G will be a three
%  dimensional array of dimensions ny x nu x M. 
%
%      
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness, 

function G = frmimo(A,B,C,s);

I = eye(size(A)); M = length(s);
ny = size(C,1); nu = size(B,2);
G = zeros(ny,nu,M);

for k=1:M
 G(:,:,k) = C*((s(k)*I-A)\B);
end;


