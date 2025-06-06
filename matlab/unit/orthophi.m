%   This generates orthonormal basis function vectors for continuous time
%   frequency domain system identification.  
%
%   Usage is
%
%   [GAMMA,NUM,DEN] = orthophi(w,poles,T.op)
%
%   where
%
%   w     =  Vector of frequencies at which measurements were conducted.
%            Element f(k) is frequency response at w(k) rad/s.
%   poles =  Vector of Poles to be used in basis functions.  
%            They are only allowed to be real valued.  Don't forget 
%	         that these should all be negative if you are specifying a stable model.
%   T     =  Sampling period in seconds.  Not used in s operator case.
%   op    =  may be 's' for Laplace, 'q' for shift or 'd' for delta operator.
%
%   GAMMA =  basis functions, one per column, evaluated at frequencies
%            specified in w.  			
%   NUM   =  If model is a linear combination of columns of GAMMA, then
%   DEN      transfer function form has a numerator which is the same linear
%            combination of the rows of NUM.  Denominator of transfer function form is DEN.
%
%
%   written by Brett Ninness, Department of EE & CE
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness

function [GAMMA,NUM,DEN] = orthophi(w,poles,T,op)

%  Now generate the regressors as specified by the poles.  This is  
%  the matrix GAMMA which also allows the frequency response of the model 
%  to be easily calculated.
     
GAMMA = []; 

%  Get frequency domain variable according to time domain operator.

if (op == 's') z = j*w;
  elseif( op == 'q') z = exp(j*w*T);
    elseif( op == 'd') z = (exp(j*w*T) - ones(size(w)))/T;  
end;  
      
%  Initialise the construction

if (op ~= 'q')
  n1 = sqrt( 2*real(-poles(1)) ); d1 = [1,-poles(1)]; NUM = [];
else
  n1 = sqrt( 1-abs(poles(1))^2 ); d1 = [1,-abs(poles(1))]; NUM = [];
end;  
  
if ( abs( imag( poles(1) ) ) > 0 )           
   error('Only real valued poles are permitted');
else
  gamma = polyval(n1,z)./polyval(d1,z);
  GAMMA = [GAMMA,gamma(:)];
  GAMMA_LAST = gamma;  num = n1;  DEN = d1;
  NUM = [NUM;zeros(1,length(poles)-length(num)),num];
end;
   
%  Now iterate through for as many poles as specified.  

for k = 2:length(poles)
  if (op ~= 'q')
    n1 = sqrt( real(-poles(k))/real(-poles(k-1)) )*[-1,-poles(k-1)'];
  else
    n1 = sqrt((1-abs(poles(k))^2)/(1-abs(poles(k-1))^2))*[poles(k-1),-1];
  end;  
  d1 = [1,-poles(k)];
  gamma = polyval(n1,z)./polyval(d1,z);
  GAMMA = [GAMMA,GAMMA_LAST(:).*gamma(:)];
  GAMMA_LAST = GAMMA_LAST(:).*gamma(:);
  num = conv(num,n1);  
  DEN = conv(DEN,d1);
  NUM = [NUM;zeros(1,length(poles)-length(num)),num];
end;  

%  Numerators above are those corresponding to different denominators.
%  Calculate numerators corresponding to common denominator.

for k=1:length(poles)-1
  x = NUM(k,:);
  for m = 1:length(poles)-k;
    x = conv(x,[1,-poles(length(poles)-m+1)]);
  end;
  NUM(k,:) = x(length(x)-length(poles)+1:length(x));
end;

   



























