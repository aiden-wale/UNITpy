%   Definition of chi squared probability density function with d degrees 
%   of freedom        
%
%	Usage y = chi(x,d);
%
%   where
%
%	x:   X axis value
%	d:   degrees of freedom
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%               		      Australia.

% Copyright (C) Brett Ninness 

function y = chi(x,d)
          
for k = 1:length(x)
 if x(k) > 0.0
  y(k) = x(k)^((d-2)/2).*exp(-0.5*x(k));
  denom = 2^(0.5*d)*gamma(0.5*d);
  y(k) = y(k)/denom;
 else
  y(k) = 0.0;         
 end;    
end;               
