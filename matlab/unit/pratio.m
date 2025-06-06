%  PRATIO: Compute the ratio of the posteriors distribution of parameters
%  given data for transfer function models.  This function is not meant
%  to be called directly by the user.  Rather, it is meant to be called
%  by POSTDIST via it employing a Metropolis-Hastings method to which
%  a handle to this function is passed.
%
%  [prat,cold] = pratio(Z,M,OPT)
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

function [prat,cold] = pratio(Z,M,OPT)

% First, set or not flag for a GOFAST, that only works for OE structures.
% It's here only for testing purposes because it is insanely faster than
% a call to VN.
 
GOFAST = 1;   % Do gofast in case we have OE structure

if GOFAST
 % If OPT.cold has a value in it, then metropolis.m is telling us that on
 % previous iteration proposal was not accepted, and hence value in
 % OPT.cold need not be recomputed
 if isempty(OPT.cold)
  num = M.thold(1:length(M.B));
  den = [1,M.thold(length(M.B)+1:end)'];
  yp = filter([zeros(1,M.delay),num],den,Z.u);   
  peold = Z.y(:)-yp(:); % Prediction error
  cold = peold'*peold/length(peold);
 else
  cold = OPT.cold;
 end;
  
 num = M.thnew(1:length(M.B));
 den = [1,M.thnew(length(M.B)+1:end)'];
 yp = filter([zeros(1,M.delay),num],den,Z.u);   
 penew = Z.y(:)-yp(:); % Prediction error
 cnew = penew'*penew/length(penew);
else   
 % This should automatically handle any model structure that VN has been
 % writtent to handle - but it will be sloooowwww.
 [cold,peold] = VN(Z,M.thold,OPT,M);   
 [cnew,penew] = VN(Z,M.thnew,OPT,M);  
end;

% Use cold, cnew, peold, penew to compute acceptance probabilities 

if OPT.rej  % Is a rejection sampling ratio requested?
 if strcmpi(OPT.dens,'gaussian') 
  prat = exp((-0.5/OPT.var)*(cnew-cold)*length(penew));   
 elseif strcmpi(OPT.dens,'uniform')    
  ptop = 1-max(abs(penew)>sqrt((3*OPT.var)));  
  pbot = (1/(sqrt(2*pi)))*exp(-0.5*M.thnew'*M.thnew);
  prat = ptop/pbot;
 end;
else % Case of a Metropolis ratio being requested
 if strcmpi(OPT.dens,'gaussian') 
  %prat = exp((-0.5/OPT.var)*(cnew-cold)*length(penew));   
  
  prat = ((sqrt(OPT.varold/OPT.var))^length(penew))*exp( -0.5*length(penew)*(cnew/OPT.var-cold/OPT.varold) );   
  
 elseif strcmpi(OPT.dens,'anneal')
  prat = exp(-(cnew-cold)/Temp);       
  if ~mod(k,10) Temp = dfac*Temp; end;
 elseif strcmpi(OPT.dens,'uniform')    
  
  knew = sqrt((3*OPT.var));
  kold = sqrt((3*OPT.varold));

  p1 = max(abs(penew)>knew);  
  prat = ((kold/knew)^length(penew))*(1-p1);
    
 end;
end;

