% STARTNL - function to initialise estimate of any nonlinear components
% of the model structure that are not initialised by the user.  This
% function is not meant to be called by a user - instead it is just an 
% auxiliary function that is used internally by other routines; 
% most importantly EST.m
%
% Usage is:
%
% M = startNL(Z,M)
%
% written by Brett Ninness, School of EE & CS
%                           University of Newcastle
%                           Australia.

% Copyright (C) Brett Ninness.


function M = startNL(Z,M)

if ~isfield(M,'in'),
 for k=1:M.nu,
  M.in(k).type = 'linear';
 end
elseif isempty(M.in),
 for k=1:M.nu,
  M.in(k).type = 'linear';
 end
else
 for k=1:M.nu,
  if ~isfield(M.in(k),'type'),
   M.in(k).type = 'linear';
  end
 end
end
if ~isfield(M,'out'),
 M.out.type = 'linear';
elseif ~isfield(M.out,'type'),
 M.out.type = 'linear';
end
if M.nu == 0,
 M.in = [];
end
if M.ny == 0,
 M.out = [];
end


%Loop over input nonlinearities to specify default values
for k=1:M.nu % Loop over all inputs setting defaults
 switch M.in(k).type,
  case 'linear'
   M.in(k).eta  = [];
   M.in(k).neta = 0;
   
  case 'hinge'
   % If hinging hyperplane non-linearity, then remove spurious specs for other types.
   if isfield(M.in(k),'upper') M.in(k).upper = [];   end;
   if isfield(M.in(k),'lower') M.in(k).lower = [];   end;
   if (~isfield(M.in(k),'eta') || (isfield(M.in(k),'eta') && isempty(M.in(k).eta)))
    M.in(k).eta = [0.05,1,-0.05,-1,-0.05,1];  % Tiny Deadzone default
   else
    if rem(length(M.in(k).eta),2)
     error('Initial guess in M.in.eta must contain pairs of parameters and hence must be of even length');
    end
   end
   M.in(k).neta = length(M.in(k).eta);
   
  case 'poly'
   % If polynomial non-linearity, then remove spurious specs for other types.
   if isfield(M.in(k),'upper') M.in(k).upper = []; end
   if isfield(M.in(k),'lower') M.in(k).lower = []; end
   if ( ~isfield(M.in(k),'eta') | (isfield(M.in(k),'eta') & isempty(M.in(k).eta) ) )
    M.in(k).eta = [1,zeros(1,11)]; % Linear function default
   elseif (length(M.in(k).eta)==1)  % Polynomial order the only specification?
    M.in(k).eta = [1;zeros(floor(M.in(k).eta)-1,1)];
   end
   M.in(k).neta = length(M.in(k).eta);
   
  case {'saturation','deadzone'}
   if isfield(M.in(k),'eta') M.in(k).eta = []; end;
   if ~isfield(M.in(k),'lower'), M.in(k).lower = []; end
   if ~isfield(M.in(k),'upper'), M.in(k).upper = []; end
   if [isempty(M.in(k).lower) isempty(M.in(k).upper)]
    M.in(k).upper = 0.5; M.in(k).lower = -0.5; M.in(k).neta = 2;
   elseif [~isempty(M.in(k).lower) ~isempty(M.in(k).upper)]
    %  Sanity check its specs for saturation and deadzone non-linearities
    if (M.in(k).lower>M.in(k).upper)
     error(sprintf('Must have M.in.lower < M.in.upper on input #%d',k));
    end;
    M.in(k).neta = 2;
   else
    if ~isempty(M.in(k).lower)
     M.in(k).upper = -M.in(k).lower; M.in(k).neta = 1;
    else
     M.in(k).lower = -M.in(k).upper; M.in(k).neta = 1;
    end;
   end
   
  otherwise
   error(['Error: Input non-linearity on input ',int2str(k),' must be one of linear, saturation, deadzone, hinge or poly']);
 end  % End of loop over all the inputs
end

%Now set defaults for output nonlinearity
switch M.out.type,
 case 'linear'
  M.out.eta  = [];
  M.out.neta = 0;
  
 case 'hinge'
  % If hinging hyperplane non-linearity, then remove spurious specs for other types.
  if isfield(M.out,'upper') M.out = rmfield(M.out,'upper');   end;
  if isfield(M.out,'lower') M.out = rmfield(M.out,'lower');   end;
  if ~isfield(M.out,'eta')  M.out.eta = [0,1,0,0,0,0];  end;  % Linear default
  M.out.neta = length(M.out.eta);
  
 case 'poly'
  % If polynomial model for input non-linearity, then remove spurious specs for other types.
  if isfield(M.out,'upper') M.out = rmfield(M.out,'upper'); end;
  if isfield(M.out,'lower') M.out = rmfield(M.out,'lower'); end;
  if ~isfield(M.out,'eta')
   M.out.eta = [0.1,1,0.1,0.1,0.01,0.01];
  elseif (length(M.out.eta)==1)  % Polynomial order the only specification?
   M.out.eta(k) = [1;zeros(floor(M.out.eta)-1,1)];
  end;
  M.out.neta = length(M.out.eta);
  
 case {'saturation','deadzone'}
  if isfield(M.out,'eta') M.out = rmfield(M.out,'eta'); end;
  if [~isfield(M.out,'lower') ~isfield(M.out,'upper')]
   M.out.upper = 0.1; M.out.lower = -0.1;   M.out.neta = 2;
  elseif [isfield(M.out,'lower') isfield(M.out,'upper')]
   %  Sanity check its specs for saturation and deadzone non-linearities
   if (M.out.lower>M.out.upper) error('Must have M.out.lower < M.out.upper'); end;
   M.out.neta = 2;
  else
   if isfield(M.out,'lower')
    M.out.upper = -M.out.lower; M.out.neta = 1;
   else
    M.out.lower = -M.out.upper; M.out.neta = 1;
   end;
  end
  
 otherwise
  error(['Error: Output non-linearity must be one of linear, saturation, deadzone, hinge or poly']);
end
