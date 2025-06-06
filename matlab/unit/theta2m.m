% THETA2M - function to convert from stacked parameter vector to model
% structure form.  This function is not meant to ever be
% called by a user - instead it is just an auxiliary function that is
% used internally by other routines; most importantly EST.m
%
% To be more exact, this function takes an initial specification M of a
% model structure together with further parameters theta.  The two are
% blended to form a final model structure G, theta is assumed to
% parameterise only those parts of G that cannot be estimated in closed
% form.
%
% Therefore, for a linear OE structure, M only specifies model
% orders, and theta is the parameters of all the dynamics.  However, for
% an ARX or FIR, all the linear dynamics would already be specified in M,
% and theta would only specify any non-linear dynamics
%
% This separation might seem artificial, but in fact is important for the
% case where we run an iterative procedure for finding estimates of
% non-linear parts, but get estimate of linear part at each iteration in
% closed form.
%
% SD is a structure with vector elements containing the standard
% deviations of estimated parameters.  These are formed from M.P, the
% covariance matrix of *all* estimated parameters, whether found in
% closed form or via an iterative procedure.
%
% Usage is:
%
% [G,SD] = theta2m(theta,M,fast)
%
% Where
%
% M     = Model structure definition in MATLAB structure.
% theta = Parameters in M stacked into vector.
% fast  = Flag, which if set to one, does not bother with conversion
%         between tf and ss forms. Default is fast=0;
%
% G     = parameters in theta put into model structure of form
%         defined by M.
% SD    = Standard deviations of parameter estimates contained in theta.
%         These are put in structure form, and are used by the routine
%         details.m
%
%
% Written by Brett Ninness: School of EE & CS
%            Adrian Wills:  University of Newcastle
%                           Australia.

% Copyright (C) Brett Ninness

function [G,SD] = theta2m(theta,M,fast)

if nargin<3 fast=1; end;  % If user doesn't specifically ask for fast version, assume not wanted.
if ~isfield(M,'par'), M.par='unknown'; end

G = M;  % Put all fixed bits in, then overwrite those that are specified by theta

% If theta is a matrix it means our purpose here is to add probablity
% density estimates (bad code AW)
addp = size(theta,2)>1;  TH=theta; theta = TH(:,1);

% If standard dev's are being requested, then get variances of parameters
if nargout>1
 try,
  if ~strcmpi(M.type,'ss'),  % At present, error estimates of pars of ss model not supported
   P = real(sqrt(abs(diag(M.P)))); P=P(:)';
  else
   P = inf*ones(size(theta));
  end
 catch
  P = inf*ones(size(theta));
 end
end

% Get number of inputs and outputs and states if NOT a greybox model

if ~strcmp(M.par,'grey'),
 if any(strcmp(M.type,{'ss','bilin','bilinear'})),
  nu = size(M.ss.B,2);  ny = size(M.ss.C,1);  nx = size(M.ss.A,1);
 else
  % How may poles/zeros in biggest denominator/numerator? (if not ss form)
  nu = M.nu; 
  ny = M.ny;
  np = max(M.nA); 
  nz = max(M.nB);
 end
end

% OK, now recover description of linear dynamics from estimated parameter vector
if strcmpi(M.type,'ss')  % State space structure is easy.
 
 switch M.par,% switch according to parametrization
  case {'grey'}
   G  = feval(M.t2m,M,theta);
   G.theta = theta;
   nu = size(G.ss.B,2);  
   ny = size(G.ss.C,1);  
   nx = size(G.ss.A,1);
   
  case {'unknown','full','ddlc','struct'}
   if M.T==0,
    G.ss.A = reshape(theta(1:nx*nx),nx,nx); idx = nx*nx;
    G.ss.B = reshape(theta(idx+1:idx+nx*nu),nx,nu); idx = idx+nx*nu;
    G.ss.C = reshape(theta(idx+1:idx+nx*ny),ny,nx); idx = idx+nx*ny;
    if ~isempty(M.ss.D) G.ss.D = reshape(theta(idx+1:idx+ny*nu),ny,nu); idx = idx+ny*nu; end;
    G.ss.Q = tril(ones(nx));
    G.ss.Q(find(G.ss.Q)) = theta(idx+1:idx+nx*(nx+1)/2);
    idx    = idx + nx*(nx+1)/2;
    G.ss.S = reshape(theta(idx+1:idx+nx*ny),nx,ny); idx = idx+nx*ny;
    G.ss.R = tril(ones(ny));
    G.ss.R(find(G.ss.R)) = theta(idx+1:idx+ny*(ny+1)/2);
    idx    = idx + ny*(ny+1)/2;
    
    Pi = [G.ss.Q zeros(nx,ny);G.ss.S' G.ss.R]*[G.ss.Q zeros(nx,ny);G.ss.S' G.ss.R]';
    G.ss.Q  = Pi(1:nx,1:nx);
    G.ss.S  = Pi(1:nx,nx+1:nx+ny);
    G.ss.R  = Pi(nx+1:nx+ny,nx+1:nx+ny);
    
   else
    G.ss.A = reshape(theta(1:nx*nx),nx,nx); idx = nx*nx;
    G.ss.B = reshape(theta(idx+1:idx+nx*nu),nx,nu); idx = idx+nx*nu;
    G.ss.C = reshape(theta(idx+1:idx+nx*ny),ny,nx); idx = idx+nx*ny;
    if ~isempty(M.ss.D) G.ss.D = reshape(theta(idx+1:idx+ny*nu),ny,nu); idx = idx+ny*nu; end;
    if ~isempty(M.ss.K) G.ss.K = reshape(theta(idx+1:idx+nx*ny),nx,ny); idx = idx+nx*ny; end;
    if ~isempty(M.ss.X1) G.ss.X1 = reshape(theta(idx+1:idx+nx),nx,1); end;
   end
      
  case 'blkdiag'
   G.ss.A=zeros(nx,nx);
   for i=1:length(theta),
    if i<=nx,
     if mod(i,2),
      G.ss.A(i,i)=theta(i);
      G.ss.A(i+1,i+1)=theta(i);
     else
      G.ss.A(i-1,i)=theta(i);
      G.ss.A(i,i-1)=-theta(i);
     end
    elseif i>nx & i<=2*nx,
     G.ss.B(i-nx,1)=theta(i);
    elseif i>2*nx & i<=3*nx,
     G.ss.C(1,i-2*nx)=theta(i);
    else
     G.ss.D=theta(i);
    end
   end
 end
 
 % Convert from state space to transfer function form if not doing fast version
 if ~fast G = sstotf(G); end;
 
elseif strcmpi(M.type,'bilin') || strcmpi(M.type,'bilinear')   % Bilinear model structure case
 switch M.par,%switch according to parametrization
  case {'unknown','full','ddlc'}
   G.ss.A=reshape(theta(1:nx*nx),nx,nx); idx = nx*nx;
   G.ss.B=reshape(theta(idx+1:idx+nx*nu),nx,nu); idx = idx+nx*nu;
   G.ss.C=reshape(theta(idx+1:idx+nx*ny),ny,nx); idx = idx+nx*ny;
   if ~isempty(M.ss.D) G.ss.D=reshape(theta(idx+1:idx+ny*nu),ny,nu); idx = idx+ny*nu; end;
   if ~isempty(M.ss.K) G.ss.K=reshape(theta(idx+1:idx+nx*ny),nx,ny); idx = idx+nx*ny; end;
   if ~isempty(M.ss.F), G.ss.F=reshape(theta(idx+1:idx+nx*nx*nu),nx,nx*nu); idx=idx+nx*nx*nu; end
   if ~isempty(M.ss.G), G.ss.G=reshape(theta(idx+1:idx+ny*nx*nu),ny,nx*nu); idx=idx+ny*nx*nu; end
   if ~isempty(M.ss.X1), G.ss.X1=reshape(theta(idx+1:idx+nx),nx,1); idx=idx+nx; end
   
  case 'struct'
   idx = 0;
   nAi = nnz(G.ss.Ai);
   if nAi>0,
    G.ss.A(logical(G.ss.Ai))=theta(idx+1:idx+nAi);
   end
   idx = nAi;
   
   nBi = nnz(G.ss.Bi);
   if nBi>0,
    G.ss.B(logical(G.ss.Bi))=theta(idx+1:idx+nBi);
   end
   idx = idx+nBi;
   
   nCi = nnz(G.ss.Ci);
   if nCi>0,
    G.ss.C(logical(G.ss.Ci))=theta(idx+1:idx+nCi);
   end
   idx = idx+nCi;
   
   nDi = nnz(G.ss.Di);
   if nDi>0,
    G.ss.D(logical(G.ss.Di))=theta(idx+1:idx+nDi);
   end
   idx = idx+nDi;
   
   nKi = nnz(G.ss.Ki);
   if nKi>0,
    G.ss.K(logical(G.ss.Ki))=theta(idx+1:idx+nKi);
   end
   idx = idx+nKi;
   
   nFi = nnz(G.ss.Fi);
   if nFi>0,
    G.ss.F(logical(G.ss.Fi))=theta(idx+1:idx+nFi);
   end
   idx = idx+nFi;
   
   nGi = nnz(G.ss.Gi);
   if nGi>0,
    G.ss.G(logical(G.ss.Gi))=theta(idx+1:idx+nGi);
   end
   idx = idx+nGi;
   
   nX1i = nnz(G.ss.X1i);
   if nX1i>0,
    G.ss.X1(logical(G.ss.X1i))=theta(idx+1:idx+nX1i);
   end
   idx = idx+nX1i;
 end
 
 
else  % Not ss form, must be some sort of poly form
 if strcmpi(M.type,'arx') % ARX form is different to other poly forms
  %  Extract numerator and denominator polys (A & B) from parameter vector G.th
  if (nu>0)   % Are we only looking at an AR model?
   index = 1; G.nB = M.nB;
   for r=1:nu  % One T/F per i/o model
    G.B(r,1:max(M.nB)+1)=[theta(index:index+M.nB(r))',zeros(1,max(M.nB)-M.nB(r))]; index = index+M.nB(r)+1;
   end;
  else
   G.B = []; G.nB = 0;
  end;
  if (M.op=='d')
   G.nA = M.nA;
   if (np>0)
    G.A = M.J(2:np+1)+theta(length(theta)-np+1:length(theta))'; G.A = [1,G.A];
   else
    G.A = 1;
   end;
  else
   G.A = theta(length(theta)-np+1:length(theta))'; G.A = [1,G.A]; G.nA = M.nA;   
  end
  
 elseif ~strcmpi(M.type,'fir')  %OE, BJ, ARMA, ARMAX cases
  %if [~strcmpi(M.type,'fir'), ~strcmpi(M.type,'arx')]
  index = 1;  % Where we are up to in stepping through theta parameter vector.
  
  % Extract parameterisation of numerator
  for r=1:nu
   % If only estimating a non-linearity, then scalar gain B(1) normalised
   % to B(1)=1 => no need to include in parameter vector, otherwise do so.
   if [~strcmp(M.in(r).type,'linear')  nz < 1]
    G.B=1; if nargout>1 SD.B = 0; end;
   else
    indices = index:index+M.nB(r); padding = zeros(1,max(M.nB)-M.nB(r));
    if strcmpi(M.op,'s')
     dB = length(G.B(r,:))-length(indices);
     G.B(r,dB+1:end)=[theta(indices)',padding];
    else  % M.op = 'q' or 'd' cases
     G.B(r,:)=[theta(indices)',padding];
    end; 
    % Check to see if standard deviation estimates on parameters are being requested
    if nargout>1 SD.B(r,:)= [P(index:index+M.nB(r)),zeros(1,max(M.nB)-M.nB(r))]; end;
    if addp  % Our purpose here in theta2m is really to add probability density estimates
     ix=1;
     for t=index:index+M.nB(r)
      [h.p,h.x] = hist(TH(t,:),50);           % Get initial pdf estimate as sample histogram
      hs=kde(h);                              % Refine via kernel density smoothed pdf estimate
      G.pb(r,ix).x=hs.x; G.pb(r,ix).p=hs.p;   % Load result into output model structure
      ix = ix+1;                              % Increment count on which element of B we are up to
     end;
    end;
    index=index+M.nB(r)+1;
   end;
  end;
  
  % Extract parameterisation of denominator
  if any(strcmpi(M.type,{'arma','armax'})),
   r=1;
   G.A(r,:)=[[1,theta(index:index+M.nA(r)-1)'],zeros(1,np-M.nA(r))];
   % Check to see if standard deviation estimates on parameters are being requested
   if nargout >1 SD.A(r,:)=[P(index:index+M.nA(r)-1),zeros(1,np-M.nA(r))]; end;
   if addp  % Our purpose here in theta2m is really to add probability density estimates
    ix=1;
    for t=index:index+M.nA(r)-1
     [h.p,h.x] = hist(TH(t,:),50);           % Get initial pdf estimate as sample histogram
     hs=kde(h);                              % Refine via kernel density smoothed pdf estimate
     G.pa(r,ix).x=hs.x; G.pa(r,ix).p=hs.p;   % Load result into output model structure
     ix = ix+1;                              % Increment count on which element of A we are up to
    end;
   end;
   index=index+M.nA(r);
   
  else %BJ or OE model, in which case A means something else
   for r=1:nu
    G.A(r,:)=[[1,theta(index:index+M.nA(r)-1)'],zeros(1,np-M.nA(r))];
    % Check to see if standard deviation estimates on parameters are being requested
    if nargout >1 SD.A(r,:)=[P(index:index+M.nA(r)-1),zeros(1,np-M.nA(r))]; end;
    if addp  % Our purpose here in theta2m is really to add probability density estimates
     ix=1;
     for t=index:index+M.nA(r)-1
      [h.p,h.x] = hist(TH(t,:),50);           % Get initial pdf estimate as sample histogram
      hs=kde(h);                              % Refine via kernel density smoothed pdf estimate
      G.pa(r,ix).x=hs.x; G.pa(r,ix).p=hs.p;   % Load result into output model structure
      ix = ix+1;                              % Increment count on which element of A we are up to
     end;
    end;
    index=index+M.nA(r);
    if any(strcmpi(M.type,{'ar','arma','arx','armax'})), break; end
   end;
  end
   
  % Extract Parameterisation of noise model
  G.C = [1,theta(index:index+length(M.C)-2)'];
  if nargout>1 SD.C = P(index:index+length(M.C)-2); end;
  if addp  % Our purpose here in theta2m is really to add probability density estimates
   ix=1;
   for t=index:index+length(M.C)-2
    [h.p,h.x] = hist(TH(t,:),50);           % Get initial pdf estimate as sample histogram
    hs=kde(h);                              % Refine via kernel density smoothed pdf estimate
    G.pc(ix).x=hs.x; G.pc(ix).p=hs.p;       % Load result into output model structure
    ix = ix+1;                              % Increment count on which element of C we are up to
   end;
  end;
  
  index = index+length(M.C)-1+isempty(M.C);
  
  if strcmpi(M.type,'bj')
   G.D = [1,theta(index:index+length(M.D)-2)'];
   if nargout>1 SD.D = P(index:index+length(M.D)-2); end;
   if addp  % Our purpose here in theta2m is really to add probability density estimates
    ix=1;
    for t=index:index+length(M.D)-2
     [h.p,h.x] = hist(TH(t,:),50);           % Get initial pdf estimate as sample histogram
     hs=kde(h);                              % Refine via kernel density smoothed pdf estimate
     G.pd(ix).x=hs.x; G.pd(ix).p=hs.p;       % Load result into output model structure
     ix = ix+1;                              % Increment count on which element of D we are up to
    end;
   end;
   index = index+length(M.D)-1+isempty(M.D);
  end;
 end; % End of check on FIR or ARX
end; % End of check on whether ss model structure or not

% Recover description of input non-linearity from estimated parameter vector
if isfield(M,'in'),
 for k=1:nu
  if ~strcmpi(M.in(k).type,'linear') %  Check if input non-linearity is to be incoporated
   if (strcmpi(M.type,'fir') || strcmpi(M.type,'arx'))
    index = 1;  % In FIR and ARX cases the linear part is not found iteratively
   end;
   if (strcmpi(M.in(k).type,'saturation') || strcmpi(M.in(k).type,'deadzone'))
    eta = theta(index:index+M.in(k).neta-1);
    if nargout>1 SD.in(k).eta = P(index:index+M.in(k).neta-1); end;
    index=index+M.in(k).neta;
    if M.in(k).neta>1
     G.in(k).upper = max(eta); G.in(k).lower = min(eta);
     if nargout>1 SD.in(k).upper = SD.in(k).eta(1); SD.in(k).lower = SD.in(k).eta(2); end;
    else
     G.in(k).upper = abs(eta); G.in(k).lower = -abs(eta);
     if nargout>1 SD.in(k).upper = SD.in(k).eta; SD.in(k).lower = SD.in(k).eta; end;
    end;
   elseif ( strcmpi(M.in(k).type,'hinge') || strcmpi(M.in(k).type,'poly') )
    eta = theta(index:index+length(M.in(k).eta)-1);
    G.in(k).eta = eta;
    if nargout>1 SD.in(k).eta = P(index:index+length(M.in(k).eta)-1); end;
    index = index+length(M.in(k).eta);
   end;
  end;
  G.in(k).type = M.in(k).type;
 end;  % End of loop over all inputs
end

% Recover description of output non-linearity from estimated parameter vector
if isfield(M,'out'),
 if ~strcmpi(M.out.type,'linear') %  Check if output non-linearity is to be incoporated
  if (strcmpi(M.out.type,'saturation') || strcmpi(M.out.type,'deadzone'))
   eta = theta(index:index+M.out.neta-1);
   if nargout>1 SD.out.eta = P(index:index+M.out.neta-1); end;
   index=index+M.out.neta;
   if M.out.neta>1
    G.out.upper = max(eta); G.out.lower = min(eta);
    if nargout>1 SD.out.upper = SD.out.eta(1); SD.out.lower = SD.out.eta(2); end;
   else
    G.out.upper = abs(eta); G.out.lower = -abs(eta);
    if nargout>1 SD.out.upper = SD.out.eta; SD.out.lower = SD.out.eta; end;
   end;
  elseif ( strcmpi(M.out.type,'hinge') || strcmpi(M.out.type,'poly') )
   eta = theta(index:index+length(M.out.eta)-1);
   G.out.eta = eta;
   if nargout>1 SD.out.eta = P(index:index+length(M.out.eta)-1); end;
   index = index+length(M.out.eta);
  end;
 end;
G.out.type = M.out.type;
end



%Finally, save theta into new model
G.theta = theta;





