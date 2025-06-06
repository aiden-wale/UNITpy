%  This function accepts a data structure and tries to interpret what type
%  of data it describes, then fills in blanks where necessary.
%
%  Usage is 
%
%  Z = startZ(Z);
% 
%  where:
%
%   Z:          Can represent time or frequency domain data as follows.
%
%   Z.type:     Can be 'time' or 'frequency' and will be used to interpret
%               the data in Z. Default is Z.type='time', but will be set
%               depending on Z.y below.
%
%   TIME:    
%      Z.y:     A matrix containing the output data.  The number of rows
%               is equal to the number of data samples, and the number of 
%               columns is equal to the number of outputs. This must be
%               supplied.
%
%      Z.u:     A matrix containing the input data.  The number of rows
%               is equal to the number of data samples, and the number of 
%               columns is equal to the number of input.
%
%      Z.t:     A vector containing the sample times (can be
%               non-equidistant). Default Z.t=0:Z.T:(N-1)*Z.T, where
%               N=number of samples.
%    
%      Z.T:     Positive number represting the sample interval. Z.T is zero
%               for non-equidistant data. Default is Z.T=1;
%
%      Z.D:     Positive number represting the integrated sampling time. 
%               Z.D must be no greater than Z.T. Z.D is zero
%               for non-equidistant data. Default is Z.D=Z.T;
%
%      Z.d:     A vector (of the same length as Z.t) corresponding to the 
%               integration time for sampled continuous-time signals. 
%               Z.d = min(diff(Z.t))*ones(1,N); by default.
%
%
%
%   FREQUENCY:
%      Z.y:     A complex valued 3D-matrix containing the output data.  The
%               first two dimensions are equal to the number of outputs and
%               inputs, respectively. The third dimension is equal to the  
%               number of data samples. This must be supplied.
%
%      Z.u:     A complex valued 3D-matrix containing the input data.  The
%               first two dimensions are equal to the number of outputs and
%               inputs, respectively. The third dimension is equal to the  
%               number of data samples. If not supplied, it is assumed that
%               Z.y represents Frequency Response Functions, for which 
%               Z.u(:,:,k) = I.
%
%      Z.w:     A vector containing the frequency points corresponding to
%               samples (can be non-equidistant). If Z.w(.) > 2*pi, then
%               Z.discrete will be set to zero. Z.w must be supplied.
%    
%      Z.T:     Positive number represting the sample interval. Z.T is zero
%               for continuous data. Default is Z.T=1;
%
%
%
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness.

function Z = startZ(Z)

% Do we have any data?
if nargin<1,
 error('No data structure supplied');
end

% Try for a quick return
if isfield(Z,'passed_startZ'),
 return;
end

%OK, if passed an empty Z, then fill it in
if isempty(Z),
 Z.u = 0;
 Z.y = 0;
end
 

% Determine if Z is a structure or a matrix
if isstruct(Z),
 % Do we have a Z.y? - MUST have this.
 if ~isfield(Z,'y'),
  error('Must have a Z.y field');
 end
else % Else Z is just a matrix so make it a structure
 [N,nin] = size(Z); 
 if N<nin, 
  Z=Z.'; 
  [N,nin] = size(Z); 
 end
 Zm = Z; 
 Z  = struct;
 y  = Zm(:,1);  
 if nin>1,
  if isreal(y),
   Z.y = y;
   Z.u = Zm(:,2:nin);
  else
   Z.y(1,1,:) = y;
   Z.w = Zm(:,2);
  end
 else
  if isreal(y),
   Z.y = y;
  else
   error('First column of data is complex valued, but there is no second column (frequency).');
  end
 end
end

% Is there any data in Z.y?
if prod(size(Z.y))==0,
 error('One or more dimsensions of Z.y is 0. Must have data to proceed.');
end

% Determine type of data (time or frequency)?
%Z.type
if ~isfield(Z,'type'),
 if isreal(Z.y), 
  Z.type='time'; 
    else
  Z.type='frequency'; 
 end
else
 switch Z.type,
  case 'time',
   if ~isreal(Z.y), 
                error('startZ:typeInconsistentReal', 'Z.type is inconsistent with Z.y (Z.y should be real valued).');
            end
  case 'frequency',
   if isreal(Z.y), 
    error('startZ:typeInconsistentComplex', 'Z.type is inconsistent with Z.y (Z.y should be complex valued).'); 
   end
  otherwise,
   error('The value in Z.type is not known.');
 end
end

% Now determine missing bits and other things depending on Z.type
switch Z.type,

%--------------------------------------------------------------------------
%   TIME DOMAIN DATA
%--------------------------------------------------------------------------
 case 'time',
  % Get output data dimensions and rotate if necessary
  [Ny,p]=size(Z.y); 
  if p>Ny, Z.y=Z.y.'; end
  [Ny,p]=size(Z.y);
  Z.Ny=Ny; 
  Z.ny=p;
  Z.nu=0; %Changed if necessary soon...
  
  if isfield(Z,'u'),
   % Get input data dimensions
   [Nu,m]=size(Z.u);
   if Nu*m==0,
    Z.u = [];
   else
    if m>Nu, Z.u=Z.u.'; end
    [Nu,m]=size(Z.u);
    Z.nu=m;
    % Chack that number of input samples equals number of output
    % samples
    if Nu~=Ny,
     error('Number of input samples does NOT equal number of output samples!');
    end
   end
  else,
   Z.u=[];
  end
  
  if isfield(Z,'T'),
   if Z.T>0, % Data should be discrete
    if isfield(Z,'t'),
     % Make sure t is equidistant
     if (abs(max(diff(Z.t)) - min( diff(Z.t))) > 1e-10),
      error('Z.t is not equidistant, but Z.T>0.')
     end
    else
     Z.t=0:Z.T:(Ny-1)*Z.T;
    end
   else % Data should be continuous    
    if ~isfield(Z,'t'),
     error('Time stamps Z.t must be supplied for non-equidistant data (i.e. if Z.T=0).');
    end
   end
   
  else % No Z.T field, so try and decide
   if isfield(Z,'t'),
    if max(diff(Z.t))-min(diff(Z.t)) > 1e-12,
     Z.T=0;
    else
     Z.T=mean(diff(Z.t));
    end
   else,
    Z.T=1;
    Z.t=0:Z.T:(Ny-1)*Z.T;
   end
  end
   
  if ~isfield(Z,'D'),
   Z.D = Z.T;
  end
  
  if isfield(Z,'d'),
   if max(d) > min(diff(Z.t)),
    error('An integration time in Z.d is greater than the time between samples in Z.t!');
   elseif length(Z.d)>1 & length(Z.d)<N,
    error('Must have the same number of entries in integration times Z.d as there are samples!')
   elseif length(Z.d)==1,
    Z.d = Z.d*ones(1,Z.Ny);
   end
  else
   Z.d = min(diff(Z.t))*ones(1,Z.Ny);
  end
  
  
%--------------------------------------------------------------------------
%   FREQUENCY DOMAIN DATA
%--------------------------------------------------------------------------
 case 'frequency'
  % Help user in case of SISO data - it's not natural to make a 3 dim array for a vector
  if (ndims(Z.y)<3)
   [M, N] = size(Z.y);
   if M > N, 
    Z.y = Z.y.';
   end
   [M, N] = size(Z.y);
   if(M < 2)
    y(1,1,:) = Z.y;
    Z=rmfield(Z,'y');
    Z.y=y;
   else
    error('startZ:lackOfDimensions', 'Z.y does not have enough  dimensions.');
   end
  end
        % Get output data dimensions and rotate if necessary
  [p,m,Ny]=size(Z.y); 
  if p>Ny | m>Ny, 
   error('Z.y needs to be transposed (see help startZ).');
  end
  Z.Ny = Ny;
  Z.ny = p;
  Z.nu = m;
  
  %Check for input data (default to Z.u(:,:,k) = I if not).
  if isfield(Z,'u'),
   % Help user in case of SISO data - it's not natural to make a 3
   % dim array for a vector
   if (ndims(Z.u)<3)
    [M, N] = size(Z.u);
    if M > N,
     Z.u = Z.u.';
    end
    [M, N] = size(Z.u);
    if(M < 2)
     u(1,1,:) = Z.u;
     Z=rmfield(Z,'u');
     Z.u=u;
    else
     %Assume a multi-input signal that is given as a m x N
     %matrix
     u(:,1,:) = Z.u;
     Z=rmfield(Z,'u');
     Z.u=u;
    end
   end
  else
   Z.u = zeros(Z.ny,Z.nu,Z.Ny);
   for k=1:Z.Ny,
    Z.u(:,:,k) = eye(Z.ny,Z.nu);
   end
  end
  
  if ~isfield(Z,'w'),
   error('Frequency domain data requires a Z.w field.');
  else % There is an omega
    if max(abs(Z.w))>pi, % Continuous model required?
     if isfield(Z,'T'),
      if max(abs(Z.w))>pi/Z.T+sqrt(eps),
       warning('Frequency range is greater than expected, i.e. max(abs(Z.w)) > pi/Z.T.');
      end
    else
      Z.T=0;
     end
    else % Should be discrete
     if ~isfield(Z,'T'),
      Z.T=1;
     end
    end
  end

%--------------------------------------------------------------------------
%   OTHERWISE DATA FORMAT NOT KNOWN
%--------------------------------------------------------------------------
 otherwise,
  error('Value in Z.type is not known.');
end

Z.passed_startZ=1;