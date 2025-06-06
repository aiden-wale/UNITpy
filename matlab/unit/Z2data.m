%  Z2DATA: Takes structure defining input-output data and extracts that
%  data together with dimensions of input and output data.
%
%  Usage is 
%
%  [y,u,ny,nu,N,Z] = Z2data(Z);
% 
%  where:
%
%   Z         : Can be in one of two forms.  The standard form is for it
%               to be a record with elements Z.y and Z.u, each of which
%               are matrices with number of rows equal to the number of
%               data samples, and number of columns equal (respectively)
%               to the number of outputs and the number of inputs.  On
%               the other hand, Z can be a matrix of the form Z = [y,u] 
%               where it is assumed that y is a column vector of output 
%               measurements and u is a matrix whose columns are the
%               input measurements; in this latter MISO models are 
%               being considered.  
%   y         : A matrix containing the output data.  The number of rows
%               is equal to the number of data samples, and the number of 
%               columns is equal to the number of outputs.
%   u         : A matrix containing the input data.  The number of rows
%               is equal to the number of data samples, and the number of 
%               columns is equal to the number of input.
%   ny        : The number of outputs, which is the column dimension of
%                Z.y if Z is a record, or it is equal to 1 if Z is a matrix;
%   nu        : The number of inputs, which is the column dimension of 
%                Z.u if Z is a record, or it is the column dimension of Z
%                minus 1 if Z is a matrix.
%   N         : Number of data samples contained in Z 
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%        		              Australia.

% Copyright (C) Brett Ninness, Adrian Wills

function [y,u,ny,nu,Ny,Z] = Z2data(Z)

%Make sure Z has passed through startZ
Z=startZ(Z);

%Extract things as already recorded
switch Z.type,
 case 'time',
  y  = Z.y;
  u  = Z.u;
  %if isempty(u), u=zeros(Z.Ny,1); end
  ny = Z.ny;
  nu = Z.nu;
  Ny = Z.Ny;
  
 case 'frequency',
  y  = Z.y;
  u  = Z.w;
  ny = Z.ny;
  nu = Z.nu;
  Ny = Z.Ny;

 otherwise,
 	error('Value in Z.type not known.');
end
