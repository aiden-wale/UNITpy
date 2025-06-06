%  Function that calls underlying subspace algorithms depending on data
%  type. E.g. equidistant discrete time data corresponds to the standard
%  and well known N4SID or CCA algorithms. Likewise, frequency domain data
%  will correspond to the frequency domain subspace algorithm.
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%                             Australia.
%

% Copyright (C) Brett Ninness

function G = subspace(Z,M,OPT)

% Make sure Z has been parsed
Z=startZ(Z);

if nargin<2,
 M   = [];
 OPT = [];
elseif nargin<3,
 OPT = [];
end

% Now switch according to data type
switch Z.type
 
 case 'time'
  G = sid(Z,M,OPT);
  
 case 'frequency'
  G = fsid(Z,M,OPT);
  
 otherwise
  error('Data type (Z.type) not known');
		
end