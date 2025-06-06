%  PTARGET: Compute the value of a given posterior density we could like
%  to sample from using a slice sampler.  This function is not meant
%  to be called directly by the user.  Rather, it is meant to be called
%  by POSTDIST if the user selects a slice sampler.
%
%  pvalue = ptarget(Z,M,OPT)
%  
%  where:
%
%   Z:          Input-Output data in one of two forms.  The standard form
%               is for it to be a record with elements Z.y and Z.u, each
%               of which are matrices with number of rows equal to the
%               number of data samples, and number of columns equal (respectively)
%               to the number of outputs and the number of inputs.
%
%   M:          Data structure which defines the model structure for
%               which the posterior distribution of the parameters or functions of
%               them is required.  Type "help est"  for a detailed
%               description.  
%
%  OPT:         Data structure which defines options for the computation
%               of the probability ratio.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%        		              Australia.

% Copyright (C) Brett Ninness

function pvalue = ptarget(Z,M,OPT)

% First, set or not flag for a GOFAST, that only works for OE structures.
% It's here only for testing purposes because it is insanely faster than
% a call to VN.
 
GOFAST = 1;   % Do gofast in case we have OE structure

% Get prediction error residuals at parameter value current
% sample point M.thnew
if GOFAST
 % If OPT.cold has a value in it, then metropolis.m is telling us that on
 % previous iteration proposal was not accepted, and hence value in
 % OPT.cold need not be recomputed

 num = M.thnew(1:length(M.B));
 den = [1,M.thnew(length(M.B)+1:end)'];
 yp = filter([zeros(1,M.delay),num],den,Z.u);   
 penew = Z.y(:)-yp(:); % Prediction error
else   
 % This should automatically handle any model structure that VN has been
 % writtent to handle - but it will be sloooowwww.
 [cnew,penew] = VN(Z,M.thnew,OPT,M);  
end;

% Use the prediction error to compute the value of the target density

 if strcmpi(OPT.dens,'gaussian') 
  cnew = penew'*penew/length(penew);  
  pvalue = ((sqrt(1/OPT.var))^length(penew))*exp( -0.5*length(penew)*cnew/OPT.var ); 
 
  pvalue = log(pvalue);
  
%  pvalue = -0.5*length(penew)*( log(OPT.var) + cnew/OPT.var ); 
  
 elseif strcmpi(OPT.dens,'uniform')    
  
  knew = sqrt((3*OPT.var));    % Uniform bound implied by variance value
  p1 = max(abs(penew)>knew);   % Zero if any one residual is outside uniform bounds
 
  if p1>0  % At least one residual violated uniform bound
   pvalue = -1e100;  % Close enough to minus infinity I guess
  else
   pvalue = -length(penew)*log(knew);
  end;
 
%  pvalue = ((1/knew)^length(penew))*(1-p1);
  
 end;


