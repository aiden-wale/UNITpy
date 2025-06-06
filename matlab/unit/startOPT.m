% STARTOPT - function to initialise estimation options in case user has been
% very sparse in setting them.  This function is not meant to ever be
% called by a user - instead it is just an auxiliary function that is
% used internally by other routines; most importantly EST.m
%
% Usage is: 
%
% OPT = startOPT(OPTin)
%
% Where
%
% OPTin = whatever has already been specified as options (may be empty)
%
% written by Brett Ninness, School of EE & CS
%            Adrian Wills   University of Newcastle
%             		        Australia.

% Copyright (C) Brett Ninness.

function OPT = startOPT(OPTin,Min)

if nargin<1, 
	OPTin=[]; 
	Min=[]; 
elseif nargin<2,
	Min=[];
end


%Default values (just add values here and they will go in automatically)
o.n      = 0;           % Amount of initial data to throw away
o.dsp    = 1;           % Flag - not set means silent running
o.miter  = 100;         % Maximum number of iterations for iterative searches
o.tol    = 1e-4;
o.lmax   = 30;
o.mdec   = 1e-9;
o.fast   = 0;
o.filt   = 0;
o.step   = 1;
o.smits  = 5;
o.cost   = 'trace';
o.subtol = sqrt(sqrt(eps));
o.delta  = 1;
o.adapt  = 1;
o.dir    = 'trust';
o.ngt    = 0;
o.nht    = 0;
o.gh     = 0;
o.saveit = 1;
o.pnum   = 100;
o.emit   = 100; 
if isfield(Min,'type'),
 if any(strcmpi(Min.type,{'ss','bilin','bilinear'})),
  o.alg   = 'gn';
  o.miter = 200;
 else
  o.alg   = 'gn';
 end
else
 o.alg    = 'gn';
end
o.basis  = 'polyb';
o.smeth  = 'bilin';
o.allP   = 1;
o.cmpgrd = 0;
o.sysnd  = 'forward';
o.gradnd = 'midpoint';
o.par    = 'ddlc';

%Fill in missing or empty fields of OPTin and put into OPT
fn=fieldnames(o);
if isempty(OPTin),
 OPT=o;
else
  OPT = OPTin;  % Preserve what is already specified.
  for i=1:length(fn),
   if ~isfield(OPT,fn{i}),
    OPT=setfield(OPT,fn{i},getfield(o,fn{i}));
   elseif isempty(getfield(OPT,fn{i}))
    OPT=setfield(OPT,fn{i},getfield(o,fn{i}));
   end
  end
end;

OPT=orderfields(OPT);



