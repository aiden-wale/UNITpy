% M2F - function to compute the frequency response of a linear time
% invariant model.  This function is not ever meant to be
% called by a user - instead it is just an auxiliary function that is
% used internally by other routines; most importantly EST.m
%
% To be more explicit, this function takes a specification M of a
% model structure and computes the input-output frequency response of the
% linear part of the model, and then adds that information to the
% structure defining the model.
%
% Usage is:
%
% G = m2f(M)
%
% Where
%
% M     = Model structure definition
%
% G     = Model structure M as above, but augmented to also include an
%         element G.G which specified the frequency response of the
%         linear part of M
%
% written by Brett Ninness, School of EE & CS
%                           University of Newcastle
%      		                   Australia.


% Copyright (C) Brett Ninness.

function G = m2f(M);

% Check if frequency domain data was passed in - allows a quick return
freqdata = 0;
if isa(M,'numeric')
 M = startZ(M);
end
if isfield(M,'type'),
 if strcmpi(M.type,'frequency'),
  freqdata = 1;
 end
elseif isfield(M,'y'),
 if ~isreal(M.y),
  freqdata = 1;
 end
end
if freqdata,
 G   = M;
 G.G = G.y;
 return;
end

% Check if time domain data was passed in
timedata = 0;
if isfield(M,'type'),
 if strcmpi(M.type,'time'),
  timedata = 1;
 end
elseif isfield(M,'y'),
 if isreal(M.y),
  timedata = 1;
 end
elseif isreal(M)
 timedata=1;
end
if timedata,
 G   = nonpar(M);
 return;
end

% Model has not come from estimation routine => populate fields
if ~isfield(M,'finishM') M = startM(M); end;

% Pass all input information through to output
G = M;

%If there is no type, then check to see if there is a G.G
if ~isfield(G,'type'),
 return;
end

% Figure out how many inputs and outputs there are
if strcmpi(M.type,'nonpar')
    return;
else
 if isfield(M,'nu'),
  nu = M.nu;
 elseif isfield(M,'B'),
  nu = size(M.B,1);
 elseif isfield(M,'ss'),
  if isfield(M.ss,'B'),
   nu = size(M.ss.B,2);
  end
 else
  nu = 1;
 end
 if isfield(M,'ny'),
  ny = M.ny;
 elseif isfield(M,'ss'),
  if isfield(M.ss,'C'),
   ny = size(M.ss.C,1);
  end
 else
  ny = 1;
 end
end

% Check for lack of input, and set to defaults
if ~isfield(G,'delay')
 G.delay = zeros(nu,1);
elseif isempty(G.delay)
 G.delay = zeros(nu,1);
end

% Set M.w and M.T to default values if not present
if ~isfield(M,'T'),
 M.T = 1;
 G.T = M.T;
end
if ~isfield(M,'w'),
 M.w = logspace(-3,log10(pi/M.T),50);
 G.w = M.w;
end

% Get appropriate discrete time frequency domain argument
if (M.op=='q'),
 ww = exp(j*M.w*M.T);
 % Compute extra phase lag implied by delays on inputs
 pdel = exp((-j*M.w(:)*M.T)*G.delay');
elseif (M.op=='d'),
 ww = (exp(j*M.w*M.T)-ones(size(M.w)))./(M.T);
 % Compute extra phase lag implied by delays on inputs
 pdel = exp((-j*M.w(:)*M.T)*G.delay');
elseif (M.op=='s'),
 ww = j*M.w;
 % Compute extra phase lag implied by delays on inputs
 pdel = exp(-j*M.w(:)*G.delay');
else
 error('M.op is not known.');
end
ww = ww(:);

G.G = zeros(ny,nu,length(M.w));
G.H = zeros(ny,nu,length(M.w));

% Set G.type to unknown if it is
if ~isfield(G,'type')
 G.type = 'unknown';
end

% Handle the case of ar, arma, arx, and armax models by setting D = A
if any(strcmpi(M.type,{'ar','arma','arx','armax'})), G.D = G.A; end

% Now handle different model types properly
switch G.type, 
 case 'ss',
  if isfield(G.ss,'K'),
   if ~isempty(G.ss.K)  % Check we have a noise model
    if G.nu>0,
     G.G=frmimo(G.ss.A,G.ss.B,G.ss.C,ww);
     if numel(G.ss.D)>0,
      G.G = G.G + G.ss.D(:,:,ones(1,length(ww)));
     end
    else
     G.G=[];
    end
    DK=eye(ny); G.H=frmimo(G.ss.A,G.ss.K,G.ss.C,ww) + DK(:,:,ones(1,length(ww)));
   else
    G.G=frmimo(G.ss.A,G.ss.B,G.ss.C,ww);
    if numel(G.ss.D)>0,
     G.G = G.G + G.ss.D(:,:,ones(1,length(ww)));
    end
   end
  else   
   G.G=frmimo(G.ss.A,G.ss.B,G.ss.C,ww);
   if numel(G.ss.D)>0,
    G.G = G.G + G.ss.D(:,:,ones(1,length(ww)));
   end
  end

  for k=1:nu,
   for m=1:ny,
    pp(1,1,:)  = pdel(:,k);
    G.G(m,k,:) = G.G(m,k,:).*pp;
   end
  end
  
  
 otherwise % Not an ss model, must be poly
  
  for k=1:nu
   if any(strcmpi(M.type,{'bj','oe'})), nua = k; else nua = 1; end
   for m=1:ny
    % Get transfer function from k'th input to m'th output
    A=G.A(nua,:,m); B=G.B(k,:,m);
    switch G.op,
     case {'q','d'}
      B = fliplr(B);
      A = fliplr(A);
      G.G(m,k,:) = polyval(B,1./ww)./polyval(A,1./ww);
      pp(1,1,:)  = pdel(:,k);
      G.G(m,k,:) = G.G(m,k,:).*pp; % Take time delay on k'th input into account
     case 's'
      G.G(m,k,:) = polyval(B,ww)./polyval(A,ww);
      pp(1,1,:)  = pdel(:,k);
      G.G(m,k,:) = G.G(m,k,:).*pp; % Take time delay on k'th input into account      
    end
   end
  end
  if isfield(G,'C') & isfield(G,'D'),
   if [~isempty(G.C) ~isempty(G.D)]  % Check we have a noise model
    for k=1:ny
     for m=1:ny
      switch G.op,
       case {'q','d'},
        % Get transfer function from k'th input to m'th output
        D=G.D(k,:,m); C=G.C(k,:,m); C=[zeros(1,length(D)-length(C)) C]; C=fliplr(C); D=fliplr(D);
        G.H(m,k,:) = polyval(C,1./ww)./polyval(D,1./ww);
        
       case 's'
        G.H(m,k,:) = polyval(G.C(k,:,m),1./ww)./polyval(G.D(k,:,m),1./ww);
      end
     end
    end
   end
  end
end

% For time series case, noise spectral factor masquerades as dynamic freq resp.
if nu<1 G.G = G.H; end;
