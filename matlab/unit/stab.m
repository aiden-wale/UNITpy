% STAB - Function which stabilises a polynomial with respect to the
% stability boundary via reflection.
%
% This function is not meant to ever be called by a user - instead it is
% just an auxiliary function that is used internally by other routines; most
% importantly EST.m
%
% Usage is: 
%
% P = stab(P,op,T)
%
% Where
%
% P  = polynomial to be stabilised 
% op = Operator that P is a polynomial in. Set to 'q' for shift and 'd' 
%      for delta.  
% T  = Sampling period in seconds.  Only used for the delta operator case
%      and default value is 1.
%
% written by Brett Ninness: School of EE & CS
%            Adrian Wills:  University of Newcastle
%         		                Australia.

% Copyright (C) Brett Ninness

function z = stab(a,op,delta);

if nargin<3 delta=1; end;  % Default sampling period is 1s if not spec'd otherwise

if op=='q'
 a=a/a(1); rts=roots(a); idx=find(abs(rts)>1); rts(idx)=1./rts(idx); z=poly(rts);
else
 bit = 1e-3;  % Roots are stable if within a circle of radius (1-bit)/delta;
 r = a(1);  a = a/r;  % Make the polynomial monic
 zz = roots(a);       % Find roots
 zz = zz + ones(size(zz))/delta; % Make the roots sit around the origin.
 R = ((1-bit)/delta)*ones(size(zz)); %  make a circle of radius of (1-bit)/delta;
 zzm=abs(zz); kk = (zzm>R); %  Which roots are outside the radius?
 zzp=angle(zz);
 zzm(kk) = R(kk); %  roots outside the radius of (1-bit)/delta get moved
 z = real(poly(zzm.*exp(j*zzp) - ones(size(zz))/delta));
end;