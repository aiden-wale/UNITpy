% M2THETA - function to convert from model structure definition to stacked
% parameter vector theta.  This function is not meant to ever be
% called by a user - instead it is just an auxiliary function that is
% used internally by other routines; most importantly EST.m
%
% Note - a key point is that it extracts only those parameters from the
% model structure that need to be found by iterative search.  Those that
% may be found in closed form (eg, parameters of linear FIR and ARX
% dynamics) are ignored.
%
% Usage is:
%
% theta = m2theta(M)
%
% Where
%
% M     = Model structure definition in MATLAB structure.
% theta = parameters in M stacked into vector.
%
% written by Brett Ninness, School of EE & CS
%            Adrian Wills   University of Newcastle
%      		                Australia.

% Copyright (C) Brett Ninness

function theta = m2theta(M);

if ~isfield(M,'par'), M.par='unknown'; end

th = [];  % Nuthin in parameter vector to start with

if strcmp(lower(M.type),'ss')  % State space structure is easy.
 switch M.par,
  case {'grey'}
   th = M.theta;
   
  case {'full','unknown','ddlc','struct'}
   if M.T==0,
    nx = size(M.ss.A,1); ny = size(M.ss.C,1);
    Qi = find(tril(ones(nx)));
    Ri = find(tril(ones(ny)));
    th = [th;M.ss.A(:);M.ss.B(:);M.ss.C(:);M.ss.D(:);M.ss.Q(Qi);M.ss.S(:);M.ss.R(Ri)];
   else
    th = [th;M.ss.A(:);M.ss.B(:);M.ss.C(:);M.ss.D(:);M.ss.K(:);M.ss.X1(:)];
   end

  case 'blkdiag'
   n   = size(M.ss.A,1); 
   for i=1:2:n, 
    th(i)=M.ss.A(i,i); 
    th(i+1)=M.ss.A(i,i+1); 
   end
   th=[th(:);M.ss.B(:);M.ss.C(:);M.ss.D(:)];
 end
 [np,nz]  = size(M.ss.A); % np=max #poles, nz= max #zeros
 [dum,nu] = size(M.ss.B); % nu=#inputs,
   
elseif strcmp(lower(M.type),'bilin') | strcmp(lower(M.type),'bilinear')
 switch M.par,
  case {'full','unknown','ddlc'}
   th=[th;M.ss.A(:);M.ss.B(:);M.ss.C(:);M.ss.D(:);M.ss.K(:);M.ss.F(:);M.ss.G(:);M.ss.X1(:)];

  case 'struct'
   th = [
         M.ss.A(logical(M.ss.Ai));...
         M.ss.B(logical(M.ss.Bi));...
         M.ss.C(logical(M.ss.Ci));...
         M.ss.D(logical(M.ss.Di));...
         M.ss.K(logical(M.ss.Ki));...
         M.ss.F(logical(M.ss.Fi));...
         M.ss.G(logical(M.ss.Gi));...
         M.ss.X1(logical(M.ss.X1i))];
 end
 
 [np,nz]  = size(M.ss.A); % np=max #poles, nz= max #zeros
 [dum,nu] = size(M.ss.B); % nu=#inputs,
else
 nu = M.nu;  % nu=#inputs, np=max #poles+1
 np = M.nA;
 nz = M.nB;
 % First put description of linear dynamics into parameter vector
 % if [~strcmp(lower(M.type),'fir'), ~strcmp(lower(M.type),'arx')]
 % Stack in parameterisation of numerators
 for r=1:nu,
  % If only estimating a non-linearity, then scalar gain B(1) normalised
  % to B(1)=1 => no need to include in parameter vector, otherwise do so.
  if (strcmp(M.in(r).type,'linear') | nz>0),
   th = [th;M.B(r,1:M.nB(r)+1)'];
  end
 end

 if any(strcmpi(M.type,{'ar','arma','arx','armax'})), nud = 1; else, nud = nu; end
 
 % Stack in parameterisation of denominators
 for r=1:nud th = [th;M.A(r,2:M.nA(r)+1)']; end;

 % Stack in noise model
 if any(strcmpi(M.type,{'arma','armax','bj'})), th = [th;M.C(2:length(M.C))']; end
 if strcmpi(M.type,'bj') th = [th;M.D(2:length(M.D))']; end

 % end;  % End of test on whether there are parameters in linear model to go into theta
end; % End of test on whether ss model structure or not

%  Put description of input non-linearity into parameter vector
nonlin = [];
for k=1:nu % Loop over all inputs
 if ~strcmp(lower(M.in(k).type),'linear') %  First Check if input non-linearity is in model
  if (M.in(k).neta == 1)  % Is a symmetric deadzone or saturation the case?
   nonlin = [nonlin;M.in(k).upper];
  else
   if isfield(M.in(k),'lower') nonlin = [nonlin;M.in(k).lower];  end;
   if isfield(M.in(k),'upper') nonlin = [nonlin;M.in(k).upper];  end;
   if isfield(M.in(k),'eta')   nonlin = [nonlin;M.in(k).eta(:)]; end;
  end;
 end;
end;

%  Put description of output non-linearity into parameter vector
if ~strcmp(lower(M.out.type),'linear') %  First Check if input non-linearity is in model
 if (M.out.neta == 1)  % Is a symmetric deadzone or saturation the case?
  nonlin = [nonlin;M.out.upper];
 else
  if isfield(M.out,'lower') nonlin = [nonlin;M.out.lower];  end;
  if isfield(M.out,'upper') nonlin = [nonlin;M.out.upper];  end;
  if isfield(M.out,'eta')   nonlin = [nonlin;M.out.eta(:)]; end;
 end;
end;

% Stack parameterisation of static nonlinearities
theta = [th;nonlin];




