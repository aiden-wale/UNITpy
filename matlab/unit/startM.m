% STARTM - function to initialise a model structure in case user has been
% very sparse in defining it.  
%
% This function tries to determine the modle type, the number of inputs and
% outputs, and the order of the various polynomials or state-space matrices
% for the system.
%
% This function is not meant to ever be
% called directly by a user - instead it is just an auxiliary function that is
% used internally by other routines; most importantly EST.m
%
% Usage is: 
%
% M = startM(Z,Min)
%
% Where
%
% Z     : Data for the estimation experiment.  Refer to est.m for
%         specification of its format. 
% Min   : Whatever has already been specified in the model structure
%
% Written by Brett Ninness  School of EE & CS
%            Adrian Wills   University of Newcastle
%                           Australia

% Copyright (C) Brett Ninness.

function M = startM(varargin)

%-------------------------- HOW IT WORKS ----------------------------------
% 1. First up, this function tries to determine if both a data set Z and an
%    initial model M are supplied.
% 2. Once this is dertermined a default model "m" is specified. This should
%    contain almost all possible variables that can be set during a call to
%    any est() routine.
% 3. Then the hard bit, determine what type of model is being asked for. If
%    M.type is specified then it will determine all further decisions. If,
%    on the other hand, it is not specified then some default decisions are
%    made based on what else has been given by the data set Z and initial
%    model M. If nothing is available, then the default model is returned.
% 4. After the model type is finalised, default values for that type are
%    populated on a case-by-case basis (via a switch).
% 5. Finally, any remaining defaults are populated and the new model
%    structure is returned.
%--------------------------------------------------------------------------

%Get the number of input arguments
numargs = length(varargin);

%------------------------------ Get Z and M -------------------------------
%Based on the number and type of input arguments, we will make certain
%default decisions
switch nargin
 case 0
  % No model or data supplied so set to defaults
  ny = 1;
  nu = 1;
  Z  = struct();
  M  = struct();
  
 case 1  % Either just Z or just M
  
  Z  = struct();
  nu = 1;
  ny = 1;
  
  % Determine M based on what was passed in
  if isa(varargin{1},'char'),
   M.type = varargin{1};
  elseif isa(varargin{1},'struct'),
   if isfield(varargin{1},'type'),
    if strcmpi(varargin{1}.type,{'time','frequency'}),
     M = [];
    else
     M  = varargin{1};
     if isfield(M,'nu'),
      nu = M.nu;
     end
     if isfield(M,'ny'),
      ny = M.ny;
     end
    end
   else
    M  = varargin{1};
    if isfield(M,'nu'),
     nu = M.nu;
    end
    if isfield(M,'ny'),
     ny = M.ny;
    end
   end
  elseif isa(varargin{1},'numeric'),
   M.nA = varargin{1};
  end
  
 otherwise  % Both Z and M passed in
  
  % Not sure what the user is asking for, do something!
  Z  = startZ(varargin{1});
  nu = Z.nu;
  ny = Z.ny;
  
  % Determine M based on what was passed in
  if isa(varargin{2},'char'),
   M.type = varargin{2};
  elseif isa(varargin{2},'struct'),
   M  = varargin{2};
  elseif isa(varargin{2},'numeric'),
   if varargin{2}>=0,
    M.nA = varargin{2};
   else
    M = struct();
   end
  end
  if isfield(M,'nu'),
   nu = M.nu;
  end
  if isfield(M,'ny'),
   ny = M.ny;
  end
end
%--------------------------------------------------------------------------


%------------------THE DEFAULT MODEL STRUCTURE-----------------------------
% All default values are here
% The default model is called "m"
% The default model looks like this........

gord    = 5;                 % Default order of G dynamics
hord    = 2;                 % Default order of H dynamics
m.type  = 'arx';             % Default type of model

% Transfer function defaults
m.A     = [];
m.B     = [];
m.C     = [];
m.D     = [];
m.nA    = gord*ones(nu,1);  % Orders for transfer function models
m.nB    = gord*ones(nu,1);
m.nC    = hord;
m.nD    = hord;

% State-space defaults
m.nx    = gord;
m.estD  = 1;
m.estK  = 1;
m.estX1 = 1;
m.estF  = 1;
m.estG  = 1;
m.par   = 'ddlc';           % Parametrization of model

% Other defaults
m.theta   = [];
m.nu      = nu;
m.ny      = ny;
m.finishM = 'finishM';

% If a sample time was given by the data, then use it
if isfield(Z,'T'),
 m.T       = Z.T;
else
 m.T       = 1;
end

% Default operator depends on the type of data (time domain = 'q', freq. dom. = 's')
if isfield(Z,'type')
 switch Z.type
  case 'frequency'
   m.op = 's';
  otherwise
   m.op = 'q';
%   if isfield(Z,'t')  % If timescale specified - could be irregular sampling
%    if (max( diff(Z.t))-min(diff(Z.t) )>1e-10)  % Is it irregular?
%     m.op = 's';  % Yes, irregular samples, must use continuous time model
%    end;
%   end;
 end
else
  m.op = 'q';
end

% If the data has an omega, then use it
if isfield(Z,'w'),
 m.w = Z.w;
else
% if strcmp(m.op,'q')
%  m.w = logspace(log10(pi/m.T/1000),log10(pi/m.T),1000);  % Frequency samples in normalised form
% elseif strcmp(m.op,'s')
%  m.w = logspace(-1,2,1000);  % Good grief! Continous time and not specified by user? Could be anything - make a wild guess
% end;
 m.w = logspace(log10(pi/m.T/1000),log10(pi/m.T),1000);  % Frequency samples in normalised form
end;
m.delay   = zeros(max(nu,1),1);

%--------------------------------------------------------------------------

%----------------------WHAT IS THE MODEL TYPE?-------------------------
% If type not specified then guess it
if ~isfield(M,'type'),
 
 % Determine what the user has supplied in terms of poly's, orders and
 % state-space matrices
 isApoly = isfield(M,'A') || isfield(M,'nA');
 isBpoly = isfield(M,'B') || isfield(M,'nB');
 isCpoly = isfield(M,'C') || isfield(M,'nC');
 isDpoly = isfield(M,'D') || isfield(M,'nD');
 if isfield(M,'ss'),
  isAss = isfield(M.ss,'A');
  isBss = isfield(M.ss,'B');
  isCss = isfield(M.ss,'C');
 else
  isAss = 0;
  isBss = 0;
  isCss = 0;
 end
 
 % Based on supplied orders and/or poly's, then try and guess type
 if nu==0,
  if  isApoly && ~isBpoly && ~isCpoly && ~isDpoly, M.type='ar';    end
  if  isApoly && ~isBpoly &&  isCpoly && ~isDpoly, M.type='arma';  end
 end
 if ~isApoly &&  isBpoly && ~isCpoly && ~isDpoly, M.type='fir';   end
 if  isApoly &&  isBpoly && ~isCpoly && ~isDpoly, M.type='arx';   end
 if  isApoly &&  isBpoly &&  isCpoly && ~isDpoly, M.type='armax'; end
 if  isApoly &&  isBpoly && ~isCpoly && ~isDpoly, M.type='oe';    end
 if  isApoly &&  isBpoly &&  isCpoly &&  isDpoly, M.type='bj';    end
 %If a D order or poly is given, then must be BJ form
 if  isDpoly                                    , M.type='bj';    end
 %If missing any poly information, then look to a SS type
 if ~isApoly && ~isBpoly && ~isCpoly && ~isDpoly,
  if isAss || isBss || isCss, M.type='ss'; end
 end
 
 %If there are multiple outputs, then default to state-space model
 if ny>1,
  M.type='ss';
 end
 
 %Fill in default type if not yet specified
 if ~isfield(M,'type'),
  M.type = m.type;
 end
end
%--------------------------------------------------------------------------


%--------------------------FILL IN DEFAULTS--------------------------------
% Do this depending on the model type.  The reader might expect this to
% be done by a case/switch statement - but for cts time tf we need to
% change type to ss, and then execute that initial bit - and matlab
% case/switch does not have C-like "fall through" behaviour.

%----------------------- AR -------------------------------------------

 if sum(strcmpi(M.type,{'ar','nar'}))
  % Has the user specified an order for the A polynomial?
  if ~isfield(M,'nA'),
   if ~isfield(M,'A'),
    M.nA = gord*ones(ny);
   else
    if floor(M.A)==M.A,
     M.nA = M.A;
    else
     M.nA = length(M.A)-1;
    end
   end
  end
  
  % Make sure that the nA and A variables have only one entry
  M.nA = M.nA(1);
  if isfield(M,'A'),
   if ~isempty(M.A),
    M.A  = M.A(1,:);
   end
  end
  M.B  = 0.0;
  M.C  = 1.0;
  M.D  = 1.0;
  
  % Set default orders for B, C and D poly's
  M.nB = 0;
  M.nC = 0;
  M.nD = 0;
  
  % Make sure the number of inputs == 0
  M.nu = 0;
  M.ny = 1;
 end;  

 %----------------------------------------------------------------------
  
 %----------------------- ARMA -----------------------------------------
 if sum(strcmpi(M.type,{'arma','narma'}))
  % Has the user specified an order for the A polynomial?
  if ~isfield(M,'nA'),
   if ~isfield(M,'A'),
    M.nA = gord;
   else
    if floor(M.A)==M.A,
     M.nA = M.A;
    else
     M.nA = length(M.A)-1;
    end
   end
  end
  % Has the user specified an order for the C polynomial?
  if ~isfield(M,'nC'),
   if ~isfield(M,'C'),
    M.nC = M.nA;
   else
    if [floor(M.C)==M.C,length(M.C)<2]
     M.nC = M.C;
    else
     M.nC = length(M.C)-1;
    end
   end
  end
  
  % Make sure that the nA, nC and A,C variables have only one entry
  M.nA = M.nA(1);
  if isfield(M,'A'),
   if ~isempty(M.A),
    M.A  = M.A(1,:);
   end
  end
  M.nC = M.nC(1);
  if isfield(M,'C'),
   if ~isempty(M.C),
    M.C  = M.C(1,:);
   end
  end
  M.B  = 0.0;
  M.D  = 1.0;
  
  % Set default orders for B and D poly's
  M.nB = 0;
  M.nD = 0;
  
  % Make sure the number of inputs == 0
  M.nu = 0;
  M.ny = 1;
 end;  
 %----------------------------------------------------------------------
  

 %----------------------- FIR ------------------------------------------
 if sum(strcmpi(M.type,{'fir','nfir'}))
  %This may also include generalised FIR models for onid
  %Has the user specified an order for the B polynomial?
  if ~isfield(M,'nB'),
   if ~isfield(M,'B'),
    M.nB = gord*ones(nu,1);
    M.nA = 0*ones(nu,1);
   else
    for i=1:nu,
     if floor(M.B)==M.B,
      M.nB(i,1) = M.B(i);
      M.nA      = 0*ones(nu,1);
     else
      M.nB(i,1) = length(M.B(i,:))-1;
      M.nA      = 0*ones(nu,1);
     end
    end
   end
  end
  %If length(M.nB)~=number of inputs, then extend it so that it does
  if length(M.nB)<nu,
   df = nu-length(M.nB);
   M.nB = [M.nB(:) ; M.nB(end)*ones(1,df)];
  else
   M.nu = length(M.nB);
  end
  %If length(M.nA)~=number of inputs, then extend it so that it does
  M.nA = 0;
  if length(M.nA)~=nu,
   df = nu-length(M.nA);
   if df>0,
    M.nA = [M.nA(:) ; M.nA(end)*ones(1,df)];
   else
    M.nA = M.nA(1:end-df);
   end
  end
  M.nC = 0;
  M.nD = 0;
  if ~isfield(M,'poles'),
   M.poles = zeros(1,max(M.nB)); %Special case of FIR from onid.m
  end
  
  %Set all remaining polynomials to default values
  if ~isfield(M,'A'), M.A = ones(nu,1); end
  M.C = 1.0;
  M.D = 1.0;
 end;
 %----------------------------------------------------------------------
  
 %----------------------- ARX ------------------------------------------
 if sum(strcmpi(M.type,{'arx','narx'}))
  %Has the user specified an order for the A polynomial?
  if ~isfield(M,'nA'),
   if ~isfield(M,'A'),
    M.nA = gord;
   elseif isempty(M.A),
    M.nA = gord;
   else
    if floor(M.A)==M.A,
     M.nA = M.A;
    else
     M.nA = length(M.A)-1;
    end
   end
  end
  %Has the user specified an order for the B polynomial?
  if ~isfield(M,'nB'),
   if ~isfield(M,'B'),
    M.nB = M.nA;
   elseif isempty(M.B),
    M.nB = M.nA;
   else
    if floor(M.B)==M.B,
     M.nB = M.B;
    else
     for i=1:nu,
      M.nB = length(M.B)-1;
     end
    end
   end
  end
  %If length(M.nB)~=number of inputs, then extend it so that it does
  if length(M.nB)<nu,
   df = nu-length(M.nB);
   M.nB = [M.nB(:) ; M.nB(end)*ones(1,df)];
  else
   M.nu = length(M.nB);
  end
  
  %Make sure that the nA and A variables have only one entry
  M.nA = M.nA(1);
  if isfield(M,'A'),
   if ~isempty(M.A),
    M.A  = M.A(1,:);
   end
  end
  M.C  = 1.0;
  M.D  = 1.0;
  
  %Set default orders for C and D poly's
  M.nC = 0;
  M.nD = 0;
 end;
 %----------------------------------------------------------------------
  
 %----------------------- ARMAX ----------------------------------------
 if sum(strcmpi(M.type,{'armax','narmax'}))
  %Has the user specified an order for the A polynomial?
  if ~isfield(M,'nA'),
   if ~isfield(M,'A'),
    M.nA = gord;
   else
    if floor(M.A)==M.A,
     M.nA = M.A;
    else
     M.nA = length(M.A)-1;
    end
   end
  end
  %Has the user specified an order for the B polynomial?
  if ~isfield(M,'nB'),
   if ~isfield(M,'B'),
    M.nB = M.nA;
   else
    for i=1:nu,
     if floor(M.B)==M.B,
      M.nB(i,1) = M.B(i);
     else
      M.nB(i,1) = length(M.B(i,:))-1;
     end
    end
   end
  end
  %Has the user specified an order for the C polynomial?
  if ~isfield(M,'nC'),
   if ~isfield(M,'C'),
    M.nC = hord;
   else
    if floor(M.C)==M.C,
     M.nC = M.C;
    else
     M.nC = length(M.C)-1;
    end
   end
  end
  %If length(M.nB)~=number of inputs, then extend it so that it does
  if length(M.nB)<nu,
   df = nu-length(M.nB);
   M.nB = [M.nB(:) ; M.nB(end)*ones(1,df)];
  else
   M.nu = length(M.nB);
  end
  
  %Make sure that the nA and A variables have only one entry
  M.nA = M.nA(1);
  if isfield(M,'A'),
   if ~isempty(M.A),
    M.A  = M.A(1,:);
   end
  end
  M.D  = 1.0;
  
  %Set default order for D poly
  M.nD = 0;
 end;
 %----------------------------------------------------------------------
  
 %----------------------- OE -------------------------------------------
 if sum(strcmpi(M.type,{'oe','noe'}))
  % Has the user specified an order for the A polynomial?
  if ~isfield(M,'nA'),
   if ~isfield(M,'A'),
    M.nA = gord;
   else
    % Loop over each input to guess order
    for i=1:nu,
     if floor(M.A)==M.A,
      M.nA(i,1) = M.A(i);
     else
      M.nA(i,1) = length(M.A(i,:))-1;
     end
    end
   end
  end
  % Has the user specified an order for the B polynomial?
  if ~isfield(M,'nB'),
   if ~isfield(M,'B'),
    M.nB = M.nA;
   else
    % Loop over each input to guess order
    for i=1:nu,
     if floor(M.B)==M.B,
      M.nB(i,1) = M.B(i);
     else
      M.nB(i,1) = length(M.B(i,:))-1;
     end
    end
   end
  end
  % If length(M.nA)~=number of inputs, then extend it so that it does
  if length(M.nA)<nu,
   df = nu-length(M.nA);
   M.nA = [M.nA(:) ; M.nA(end)*ones(1,df)];
  else
   M.nu = length(M.nA);
  end
  % If length(M.nB)~=number of inputs, then extend it so that it does
  if length(M.nB)<nu,
   df = nu-length(M.nB);
   M.nB = [M.nB(:) ; M.nB(end)*ones(1,df)];
  else
   M.nu = length(M.nB);
  end
  
  % Set default orders for C and D poly's
  M.nC = 0;
  M.nD = 0;
  M.C  = 1.0;
  M.D  = 1.0;
 end; 
 %----------------------------------------------------------------------
  
 %----------------------- BJ -------------------------------------------
 if sum(strcmpi(M.type,{'bj','nbj'}))
  % Has the user specified an order for the A polynomial?
  if ~isfield(M,'nA'),
   if ~isfield(M,'A'),
    M.nA = gord;
   else
    % Loop over each input to guess order
    for i=1:size(M.A,1),
     if floor(M.A)==M.A,
      M.nA(i,1) = M.A(i);
     else
      M.nA(i,1) = length(M.A(i,:))-1;
     end
    end
   end
  end
  % Has the user specified an order for the B polynomial?
  if ~isfield(M,'nB'),
   if ~isfield(M,'B'),
    M.nB = M.nA;
   else
    % Loop over each input to guess order
    for i=1:size(M.B,1),
     if floor(M.B)==M.B,
      M.nB(i,1) = M.B(i);
     else
      M.nB(i,1) = length(M.B(i,:))-1;
     end
    end
   end
  end
  % Has the user specified an order for the C polynomial?
  if ~isfield(M,'nC'),
   if ~isfield(M,'C'),
    if isfield(M,'D') if (floor(M.D)==M.D) hord = M.D; end; end;
    M.nC = hord;
   else
    if [floor(M.C)==M.C,length(M.C)<2]
     M.nC = M.C;
    else
     M.nC = length(M.C)-1;
    end
   end
  end
  % Has the user specified an order for the D polynomial?
  if ~isfield(M,'nD'),
   if ~isfield(M,'D'),
    M.nD = M.nC;
   else
    if floor(M.D)==M.D,
     M.nD = M.D;
    else
     M.nD = length(M.D)-1;
    end
   end
  end
  % If length(M.nA)~=number of inputs, then extend it so that it does
  if length(M.nA)<nu,
   df = nu-length(M.nA);
   M.nA = [M.nA(:) ; M.nA(end)*ones(1,df)];
  else
   M.nu = length(M.nA);
  end
  % If length(M.nB)~=number of inputs, then extend it so that it does
  if length(M.nB)<nu,
   df = nu-length(M.nB);
   M.nB = [M.nB(:) ; M.nB(end)*ones(1,df)];
  else
   M.nu = length(M.nB);
  end
 end; 
 %----------------------------------------------------------------------
  
  
 %----------------------- SS -------------------------------------------
 if sum(strcmpi(M.type,{'ss','nss','bilin','bilinear'}))
  % Determine state order nx
  if ~isfield(M,'nx'),
   if isfield(M,'ss'),
    if isfield(M.ss,'A')
     M.nx = size(M.ss.A,1);
    end
   end
  end
  
  if ~isfield(M,'nx'),
   if isfield(M,'nA'),
    M.nx = M.nA;
   elseif isfield(M,'A'),
    if floor(M.A)==M.A,
     M.nx = max(max(M.A));
    else
     M.nx = m.nx;
    end
   else
    M.nx = m.nx;
   end
  end
  
  % if a system has been given, then make sure nx, nu, and ny are set
  % correctly
  if isfield(M,'ss'),
   if isfield(M.ss,'A'),
    if ~isempty(M.ss.A),
     if M.nx~=size(M.ss.A,1),
      M.nx = size(M.ss.A,1);
     end
    end
   end
   if isfield(M.ss,'B'),
    if ~isempty(M.ss.B),
     if nu~=size(M.ss.B,2),
      M.nu = size(M.ss.B,2);
     end
    end
   else
    M.nu = 0;
   end
   if isfield(M.ss,'C'),
    if ~isempty(M.ss.C),
     if ny~=size(M.ss.C,1),
      M.ny = size(M.ss.C,1);
     end
    end
   else
    M.ny = 0;
   end
  end
    
  % Make sure that if a continuous-time model is asked for then we set
  % M.par = 'full';
  discrete = 1;
  if isfield(M,'op'),
   if strcmpi(M.op,'s'),
    discrete = 0;
   end
  end
  if ~isfield(M,'par'),
   if ~discrete,
    M.par = 'full';
   end
  end
 end;       
 %----------------------------------------------------------------------
 %
  
 %----------------------- NONPAR ---------------------------------------
 if strcmpi(M.type,'nonpar') 
 end;
  
 %----------------------------------------------------------------------
  
 %----------------------- STATIC ---------------------------------------
 if strcmpi(M.type,'static')
  M.nA = 0;
  M.nB = 0;
  M.nC = 0;
  M.nD = 0;
  M.A  = 1.0;
  M.B  = 0.0;
  M.C  = 1.0;
  M.D  = 1.0;
 end;
 %----------------------------------------------------------------------
  
 %----------------------- NONLINEAR STATE-SPACE ------------------------
 if strcmpi(M.type,'nlss')
   M.op = 'q';
   M.finishM = 'nlssfinish';
 end;       
        
%  %----------------------- DEFAULT --------------------------------------
% otherwise
%  M.type = m.type;

%------------------Fill in bits of M that come from data-------------------
if isfield(Z,'type'),
 switch Z.type,
  case 'frequency',
   if ~isfield(M,'w'),
    M.w = Z.w;
   end
   if isfield(M,'T')
    if M.T*max(M.w)>pi+1000*eps,
     warning('M.T and Z.w are not compatible for M.op="q", resetting M.T = pi/max(M.w);')
     M.T = pi/max(M.w);
    end
   else
    M.T = pi/max(M.w);
   end
   
  case 'time'
   if ~isfield(M,'T'),
    M.T = Z.T;
   end
   if ~isfield(M,'w')
    if isfield(Z,'w')
     M.w = Z.w;
    else
     M.w = logspace(log10(pi/M.T)-4,log10(pi/M.T),2000);
    end
   end
 end
end
%--------------------------------------------------------------------------

%------------------Fill in missing or empty fields of M--------------------

if ~strcmpi(M.type,'nlss')  %If M.typeis nonlinear state-space then don't fill in defaults
 fn=fieldnames(m);
 if isempty(M),
  M=m;
 else
  for i=1:length(fn),
   if ~isfield(M,fn{i}),
    M = setfield(M,fn{i},getfield(m,fn{i}));
   elseif isempty(getfield(M,fn{i}))
    M = setfield(M,fn{i},getfield(m,fn{i}));
   end
  end
 end
    
 %Make sure the delays are of the correct orientation and size.
 M.delay = M.delay(:);
 if size(M.delay,1) ~= M.nu,
  M.delay = [M.delay;zeros(M.nu-size(M.delay,1),1)];
 end
  
 %Initialise Hammerstein and Wiener NL blocks
 M = startNL(Z,M);

end
    
M = orderfields(M);


