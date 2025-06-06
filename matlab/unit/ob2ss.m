%   This function converts a model parameterised in orthonormal form as
%
%   G(z) = x_1*B_1(z) + x_2*B_2(z) + ... + x_n*B_n(z)
%   
%   where the B_k(z) are the general Kautz orthonormal bases with poles at
%   p_1, p_2,....,p_n into a state space form [A,B,C,0].
%
%   Input variables expected are:
%   
%   theta:  theta = [x_1,x_2,....,x_n] is the vector of parameters 
%           defining the model with respect to the orthonormal basis.
%
%   poles:  poles = [p_1,p_2,...,p_n] are the poles used to define the Kautz
%           orthonormal basis function parameterising the model.
%
%   Output variables provided are:
% 
%   A,B,C:  State space realisation such that C(zI-A)^(-1)B = G(z).
%
%   This routine is meant to be used in conjuction with the onid estimation
%   routine.  For example
%    
%   >> [THETA,PHI,GAMMA] = onid([y(:),u(:)],poles,w,1);
%   >> [A,B,C] = ob2ss(THETA,poles);
%   
%   will give an estimated model in state space form.   
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness

function [A,B,C] = ob2ss(theta,poles)

p = length(poles);
A = diag(poles);
B = [1;zeros(p-1,1)];
C = eye(p,p);
pp = -poles;
wun = ones(size(poles));
L = sqrt(diag(wun - abs(poles).^2));

%  Run down first sub-diagonals putting in the correct entries.
for m=2:p
  for n = 1:m-1
    A(m,n) = 1-poles(m)*poles(m-1);
    C(m,n) = prod(pp(n:m-1));
  end;
end;

%  For the case of the A matrix re-do 3rd and lower sub-diagonals putting in
%  some product terms.
for m=3:p
  for n = 1:m-1
    A(m,n) = A(m,n)*prod(pp(n:m-2));
  end;
end;

C = theta(:)'*L*C;




