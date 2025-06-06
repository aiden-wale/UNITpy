% UNTREND - Function removes simple trends from data
%
% X = untrend(Z,ord,wc)
%
% Where
%
% Z   = Data to be detrended
% ord = Polynomial order of detrending - default is ord=0 which implies
%       that only constant offsets are removed. If specified as the
%       string 'hp', then UNTREND tries to remove DC offsets only by high
%       pass filtering, with optional cutoff specified by third argument.
% wc  = High pass filter cut-off frequency *normalised* for unity
%       sampling period.   Default is wc=0.01*pi;
% X   = Detrended data
%
% written by Brett Ninness, Department of EE & CE
%                           University of Newcastle
%                           Australia.

% Copyright (C) Brett Ninness.

function X = untrend(Z,ord,wc);

if nargin<2 ord = 0; end;

% Extract sizes of input and output from data matrix
[y,u,ny,nu,Ny] = Z2data(Z);

if strcmp(lower(ord),'hp')  % Has user specified hp filter?
 if nargin<3 wc = 0.01*pi; end;  % Default HP filter cutoff
% b = 4*conv([1,-1],[1,-1]);  % Specify 2nd order Butterworth
% a2=wc^2+2*sqrt(2)*wc+4; a1=2*wc^2-8; a0=wc^2-2*sqrt(2)*wc+4;
% a = [a2,a1,a0]/a2; b = b/a2;
 [b,a]= butter(2,wc/pi,'high');
 X.y = filter(b,a,y); X.u = filter(b,a,u);
else % Otherwise remove polynomial trend of given order
 % Detrend all the outputs
 idx = 0:1:Ny-1; idx=idx(:);
 for k=1:ny
  % Compute regressors
  ph = [ones(Ny,1),(idx*ones(1,ord)).^(ones(Ny,1)*(1:1:ord))];
  th = ph\y(:,k);
  X.y(:,k) = y(:,k)-ph*th;
 end;
 
 % Now do all the inputs
 for k=1:nu
  % Compute regressors
  ph = [ones(Ny,1),(idx*ones(1,ord)).^(ones(Ny,1)*(1:1:ord))];
  th = ph\u(:,k);
  X.u(:,k) = u(:,k)-ph*th;
 end;
end;