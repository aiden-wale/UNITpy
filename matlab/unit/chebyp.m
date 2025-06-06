%   Auxiliary function, used by frequency domain estimation algorithms,
%   and which calculates modified Chebychev polynomials that may be used
%   for parameterising numerator and denominators.  They are modified by
%   flipping signs of alternate co-efficients so that with complex valued
%   arguments (normalised to have magnitude < 1) magnitude of polynomials
%   is less than one.  All this is due to Bayards work in Automatica 
%
%   Usage is  X = chebyp(n)
%
%   where
%
%   n = maximum order of modified Chebychev polynomial
%   X = matrix with rows equal to co-effs of modified 
%       Chebychev polynomials, highest order to lowest order
%       do-efficients ordered left to right.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%        		              Australia.

% Copyright (C) Brett Ninness.
 
function X = chebyp(n)

X = zeros(n,n); X(1,1) = 1; if (n>1) X(2,2) = 1; end;
%  Use Recursions to generate poly coeffs as rows of matrix;
for k=3:n; X(:,k) = [0;2*X(1:n-1,k-1)] - X(:,k-2);  end; X = fliplr(X');  
%  Flip some signs to make them vary properly for complex data application
Q = ones(size(X)); m = 2; 
while (m<=n-1)
  if (2*floor(m/2) == m) k=2; else k=3; end;
  while (k<=n-1) Q(m+1,n-k) = -1; k = k+4; end;
  m = m+1;
end;
X = X.*Q;
