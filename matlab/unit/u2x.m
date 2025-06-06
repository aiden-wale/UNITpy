% U2X - function to pass an input u through a specified non-linearity X so
% as to generate a new signal x = X(u).  Also, this function will generate
% the signal z which is the sensitivity (derivative) of X(u) with respect to
% its paramaterisation.  This function is not meant to ever be called by a
% user - instead it is just an auxiliary function that is used internally by
% other routines; most importantly EST.m
%
% Usage is:
%
% [x,z,w] = u2x(u,M)
%
% Where
%
% u     = vector of input signal.
% M     = Definition of model structure from which the parameterisation
%         of X is inferred.
% x     = Signal x = X(u), with each column of u being associated
%         conformally with the columns of u
% z     = Matrix, with each column being a derivative of X(u) with
%         respect to the parameters defining X(u).  This matrix is arranged
%         in blocks, with the k'th block corresponding to the the k'th
%         column u(:,k) of u, and the block containing M.in(k).neta columns.
% w     = Vector which is the derivative of X(u) with respect u, and with
%         each column being associated conformally with the columns of u.
%
% Written by Brett Ninness, School of EE & CS
%                           University of Newcastle
%                           Australia.

% Copyright (C) Brett Ninness.

function [x,z,w] = u2x(u,M)

% Figure out dimensions of input signals
[m,nu] = size(u); if (m<nu) u=u'; [m,nu] = size(u); end;
if (nargout>1) div = 1; else div = 0; z = []; end;

% z is derivative of input non-linearity with respect to its parameters
sumit = 0; for k=1:nu sumit = sumit+M.in(k).neta; end;
z = zeros(length(u(:,k)),sumit);
zind = 1;  % Where we are up to in building up matrix of derivatives

for k=1:nu
 if strcmpi(M.in(k).type,'saturation')
  x(:,k)=sat(u(:,k),M.in(k).lower,M.in(k).upper,1);
  if div  % Only do this if derivatives are asked for
   if (M.in(k).neta==2)  % Is saturation non-symmetric?
    z(:,zind) = z(:,zind) + (u(:,k) - M.in(k).lower*ones(size(u(:,k)))<0).*(M.in(k).upper>M.in(k).lower);
    z(:,zind+1) = z(:,zind+1) + (u(:,k) - M.in(k).upper*ones(size(u(:,k)))>0).*(M.in(k).upper>M.in(k).lower);
    w(:,k) = ((u(:,k) > M.in(k).lower) & (u(:,k) < M.in(k).upper));
    zind = zind+2;
   else  % Symmetric saturation case
    eta = M.in(k).upper;
    z(:,zind) = z(:,zind) + ( (u(:,k) - eta*ones(size(u(:,k))) ) >0 ).*(u(:,k)>zeros(size(u(:,k))));
    z(:,zind) = z(:,zind) - ( (u(:,k) + eta*ones(size(u(:,k))) ) <0 ).*(u(:,k)<zeros(size(u(:,k))));
    w(:,k) = ((u(:,k) > -eta) & (u(:,k) < eta));
    zind = zind+1;
   end;
  end;
 elseif strcmpi(M.in(k).type,'deadzone')
  x(:,k)=dzone(u(:,k),M.in(k).lower,M.in(k).upper);
  if div  % Only do this if derivatives are asked for
   z(:,zind:zind+M.in(k).neta-1) = -ones(length(u(:,k)),M.in(k).neta);
   if (M.in(k).neta==2)  % Is deadzone non-symmetric    
    z(:,zind)   = -double(u(:,k) < M.in(k).lower);
    z(:,zind+1) = -double(u(:,k) > M.in(k).upper);
    w(:,k) = double((u(:,k) < M.in(k).lower) | (u(:,k) > M.in(k).upper));
    zind = zind+2;
   else  % Symmetric deadzone case
    eta = M.in(k).upper;
    z(:,zind) = double(u(:,k) < -eta) - double(u(:,k) > eta);
    w(:,k) = double((u(:,k) < -eta) | (u(:,k) > eta));
    zind = zind+1;
   end;
  end;
 elseif strcmpi(M.in(k).type,'hinge')
  [x(:,k),w(:,k)]=hinge(u(:,k),M.in(k).eta);
  if div  % Only do this if derivatives are asked for
   wun = ones(size(u(:,k)));  z = zeros(length(u(:,k)),length(M.in(k).eta));
   z(:,zind) = wun; z(:,zind+1) = u(:,k);
   for r = 1:length(M.in(k).eta)/2-1
    a1 = M.in(k).eta(2*r+1); a2 = M.in(k).eta(2*r+2);
    if (abs(a2)>eps) breakpoint = -a1/a2; else breakpoint = 1e20; end;
    index = u(:,k) > breakpoint*wun;
    z(logical(index),zind+2*r) = wun(logical(index));
    z(logical(index),zind+2*r+1) = u(logical(index),k);
   end;
  end;
  zind = zind+M.in(k).neta;
 elseif strcmpi(M.in(k).type,'poly')
  [x(:,k),z(:,zind:zind+M.in(k).neta-1),w(:,k)]=polynom(u(:,k),M.in(k).eta);
 elseif strcmpi(M.in(k).type,'linear')
  x(:,k) = u(:,k); w(:,k) = ones(size(u(:,k)));
 else
  error('Specified nonlinearity is not one of saturation, deadzone, linear, hinge, poly')
 end;
end; % End of loop over all possible inputs







