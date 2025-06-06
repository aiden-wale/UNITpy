%  delta operator version of fstab.  Function stabilises a MONIC delta
%  operator polynomial with respect to the stability boundary via
%  reflection.
%  Usage is
%
%  z = delstab(a,delta);
%
%  where a is the polynomial to be stabilised, delta is the sampling period,
%  and z is the returned stabilised polynomial.
%
%  Written by Brett Ninness: School of Elec. Eng. and Comp. Sci.
%                            University of Newcastle
%                            Australia

% Copyright (C) Brett Ninness

function z = delstab(a,delta);

bit = 1e-3;   % Roots must be in circle of radius (1-bit)/delta;

% make the polynomial monic
r = a(1);  a = a/r;
zz = roots(a);

% make the roots sit around the origin
zz = zz + ones(size(zz))/delta;

%  make a circle of radius of (1-bit)/delta;
R = ((1-bit)/delta)*ones(size(zz));

%  Which roots are outside the radius?
zzm = abs(zz);  zzp = angle(zz); kk = (zzm>R);

%  roots outside the radius of (1-bit)/delta get moved
%zzm(kk) = R(kk)/delta./zzm(kk);

zzm(kk) = R(kk);
z = real(poly(zzm.*exp(j*zzp) - ones(size(zz))/delta));

