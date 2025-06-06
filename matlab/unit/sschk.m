%  SSCHK: Takes a model structure and dimensions of a multivariable data 
%  record and sanity checks them, plus also manipulates them so as to be 
%  appropriate for an underlying mex file implementation of an EM algorithm.
%  There is no conceivable reason why a user should ever call this function.
%  It is designed to be completely internal to the workings of other routines.
%
%  Usage is 
%
%  [m,nx] = sschk(M,nu,ny);
% 
%  where:
%
%   M         : An initial state space model structure specification.  
%               Only the terms M.ss.A -- M.ss.G are of interest.
%   m         : An output model structure that is for purely internal 
%               use with regard to underlying mex files.  Basically, it is 
%               the M.ss.A--M.ss.G terms in M with the .ss. stripped out, 
%               also with some sanity checking.
%   nx:       : The state dimension of the model structure specified in M & m.
%   ny        : The number of outputs, which is the column dimension of
%               Z.y if Z is a record, or it is equal to 1 if Z is a matrix;
%   nu        : The number of inputs, which is the column dimension of 
%               Z.u if Z is a record, or it is the column dimension of Z
%               minus 1 if Z is a matrix.
%
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%                     		  Australia.

% Copyright (C) Brett Ninness.
 
function [M,nx] = sschk(M,nu,ny);

% Get number of states in model
if ~isfield(M,'A') error('System matrix A must be supplied'); else nx=size(M.A,1); end;

% Check for presence and appropriate dimensions of other system matrices
if ~isfield(M,'C')
 error('Initial estimate M.ss.C must be supplied');
elseif size(M.C,1)~=ny
 error('Height of M.ss.C matrix is inconsistent with number of outputs');
end

if ~isfield(M,'B')
 if nu>0 error('Initial estimate M.ss.B must be supplied'); else M.B = zeros(nx,nu); end;
elseif [nu>0,size(M.B,2)~=nu] error('Width of M.ss.B matrix is inconsistent with number of inputs'); end;

if ~isfield(M,'D')
 if nu>0 error('Initial estimate for M.ss.D must be supplied');
  else M.D = zeros(ny,nu); end;
elseif size(M.D)~=[ny,nu]
  error('Dimension of M.ss.D matrix is inconsistent with M.ss.B and M.ss.C');
end

% Add in defaults for bilinear component specifications if required
if strcmp(M.type,{'bilin','bilinear'})
 if ~isfield(M,'F') M.F = zeros(nx,nx*nu);
  elseif isempty(M.F) M.F = zeros(nx,nx*nu);
   elseif size(M.F)~=[nx,nx*nu],
    error('The size of M.ss.F is not consistent with M.ss.A and M.ss.B');
 end;

 if ~isfield(M,'G') M.G = zeros(ny,nx*nu); 
  elseif isempty(M.G) M.G = zeros(ny,nx*nu); 
   elseif size(M.G)~=[ny,nx*nu]
    error('The size of M.ss.G is not consistent with M.ss.C, M.ss.A and M.ss.B');
 end
end;

% Initialise estimates of state and measurement noise covariance matrices (if necessary)
if ~isfield(M,'Q')  M.Q = 10*eye(nx);
 elseif isempty(M.Q) M.Q = 10*eye(nx);
 else 
 if size(M.Q)~=size(M.A) 
  error('The covariance matrix Q is not consistent with A'); 
 end;
end

if ~isfield(M,'S') M.S = zeros(nx,ny); 
 elseif isempty(M.S) M.S = zeros(nx,ny); 
  else if size(M.S)~=size(M.C') error('The S matrix is not consistent with A and C'); end;
end

if ~isfield(M,'R') M.R = 0.1*eye(ny);
 elseif isempty(M.R) M.R = 0.1*eye(ny);
  else if size(M.R)~=[ny,ny] error('The covariance matrix R is not consistent with C'); end;
end

% Initialise estimates of initial state and its covariance (if necessary)
if ~isfield(M,'P1') M.P1 = 10*eye(nx);
 elseif isempty(M.P1) M.P1 = 10*eye(nx);
  else if size(M.P1)~=size(M.A) error('The covariance matrix P1 is not consistent with A'); end;
end

if ~isfield(M,'mu') M.mu = zeros(nx,1);
 elseif isempty(M.mu),
  else if size(M.mu)~=[nx,1] error('The initial state mu is not consistent with A'); end;
end
