% Used by postdist.m to provide an alternate means to numerically compute
% integrals using Simpson's method.  There is no conceivable reason why
% user would ever call this.
%
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.
%

% Copyright (C) Brett Ninness


function p = sumd(arg,idx,count)

global P;  % Used to keep track of running sum of integral, hence must be global

if arg.shuffle(idx) <= length(arg.num)
 x = arg.xb(arg.shuffle(idx)).x;
else
 x = arg.xa(arg.shuffle(idx)-length(arg.num)).x;
end;

if idx < length(arg.den)+length(arg.num)-1
 for k=1:length(x)
   if arg.shuffle(idx)<=length(arg.num)
    arg.num(arg.shuffle(idx)) = x(k); 
   else
    arg.den(arg.shuffle(idx)-length(arg.num)+1) = x(k);
   end;
   sumd(arg,idx+1,count);  % Descend another level in recursive for loops
 end;
else
 yp = filter(arg.num,arg.den,arg.u);   % Get predictor
	pe = arg.y(:)-yp(:); % Prediction error
	if strcmp(lower(arg.dens),'gaussian') % Gaussian density comparison
  logp = -(0.5*length(pe)*log(2*pi*arg.v)) -(0.5/arg.v)*pe(:)'*pe(:); 
  p=exp(logp); if isinf(p) p = 1e5; end;
%  if [arg.den(2) > -0.799, arg.den(2) < -0.797, arg.num(1)> 0.201 , arg.num(1) < 0.203] keyboard; end;
%	 p = ( (1/sqrt(2*pi*arg.v))^length(pe) )*exp(-(0.5/arg.v)*pe(:)'*pe(:)); 
	elseif strcmp(lower(arg.dens),'uniform')
	 p= 1-max( abs(pe)>sqrt(3*arg.v) );  
	end;
 if isnan(p) p=0; warning('Warning - pint compututed NAN for an x-axis range value'); end;
 
 if arg.shuffle(idx) <= length(arg.num)
  P.pb(arg.shuffle(idx)).p(count) = P.pb(arg.shuffle(idx)).p(count) + p;
 else
  P.pa(arg.shuffle(idx)-length(arg.num)).p(count) = P.pa(arg.shuffle(idx)-length(arg.num)).p(count) + p;  
 end;
end;


