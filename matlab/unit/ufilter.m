% This function takes a numerator and denominator polynomial in either
% backward shift or inverse delta operator, and calls underlying filter
% routines based on the operator. It is not intended to be called by the
% user, but rather from other routines in the toolbox.
%
% Function use:
%
%       y = ufilter(b,a,u,M);
%  
% where:
%
%        b:  numerator polynomial
%        a:  denominator polynomial
%        u:  input time sequence
%        y:  output time sequence
%        M:  model structure containing at least M.op, which can be either 
%            'q' for backward shift operator, or, 'd' for inverse delta
%            operator, i.e. d = (q-1)/Delta and the polynomials are in 
%            d^{-1}.
%
%    Written by Adrian Wills,  Department of EE & CE
%               Brett Ninness  University of Newcastle
%                     		   Australia.

% Copyright (C) Brett Ninness

function y=ufilter(b,a,u,M,x)

if nargin<3,
 error('ufilter requires at least 3 input arguments');
end

if nargin<4,
 M.op = 'q'; 
 M.T  = 1;
 x    = [];
end

if nargin<5,
 x  = [];
end

% Call appropriate filter based on operator.

if (M.op=='d')
 y=delfilter(b,a,u,M.T);
else
 y=filter(b,a,u);
end;
