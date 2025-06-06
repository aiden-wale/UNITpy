%   This function provides a mapping from Model type (M.type) and Data type
%   (Z.type) to different initialisation and estimation algorithms within the toolbox.
%   
%   The intention is to provide one definitive place where the functionality
%   of est.m is specified. So that, adding new Model types or Data types is a
%   matter of supplying the underlying routines and then updating estmap.m to
%   tell est.m what to do with them.
%
%   This routine is never meant to be called directly by the user - it is 
%   called by est.m 
%   
%   Usage:   ep = estmap(Z,M,OPT);
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%             		          Australia.


% Copyright (C) Brett Ninness

function ep = estmap(Z,M,OPT);

% Functions will be called in this order from est.m
% ep.startM  = 'startM'; %default unless otherwise specified
ep.startG  = '';
ep.startH  = '';
ep.startNL = '';
ep.alg     = '';
ep.finishM = 'finishM';

if nargin<3,
 warning('estmap: 3 inputs should be supplied');
 return;
end

if isempty(M),
 return;
elseif ~isstruct(M),
 return;
elseif ~isfield(M,'type'),
 return;
end

% Setup strings so that we can print model equations
if isfield(M,'in'),
 nli = 0;
 for i=1:M.nu,
  if ~strcmpi(M.in(i).type,'linear'),
   nli = 1;
  end
 end
 if nli,
  ham = '';
  if M.nu>1,
   inp = 'x_i(t)';
   for i=1:M.nu,
    ham = [ham sprintf('   x_%i(t) = %s(u_%i(t))\n',i,M.in(i).type,i)];
   end
  else
   ham = sprintf('     x(t) = %s(u(t))\n',M.in(i).type);
   inp = 'x(t)';
  end
 else
  ham = '';
  if M.nu>1,
   inp = 'u_i(t)';
  else
   inp = 'u(t)';
  end
 end
end

if isfield(M,'out'),
 if ~strcmpi(M.out.type,'linear'),
  out = 'z(t)';
  wen = sprintf('     y(t) = %s(z(t))',M.out.type);
  nlo = 1;
 else
  out = 'y(t)';
  wen = '';
  nlo = 0;
 end
else
 out = 'y(t)';
 wen = '';
 nlo = 0;
end

% Switch according to model type
switch lower(M.type),
 %----------------------------------------------------------------------
 %  NONPAR
 %----------------------------------------------------------------------
 case 'nonpar',
  ep.alg='nonpar';
  ep.modelEquations = sprintf('G(%s) = U(%s) / Y(%s)',M.op,M.op,M.op);

 %----------------------------------------------------------------------
 %  ARX
 %----------------------------------------------------------------------
 case {'ar','arx','farx'}
  switch Z.type
   case 'time'
    ep.alg='barx';
    for k=1:M.nu,
     if ~strcmpi(M.in(k).type,'linear'),
      ep.startNL = 'startNL';
     end
    end
    for k=1:M.ny,
     if ~strcmpi(M.out(k).type,'linear'),
      ep.startNL = 'startNL';
     end
    end
   case 'frequency'
    ep.alg='farx';
  end

  if strcmpi(M.type,'ar'),
   cr = sprintf('\n');
   s1 = sprintf('                          ');
   s2 = sprintf('     A(%s)y(t) = e(t)      ',M.op);
   s3 = sprintf('                          ');
   s4 = ['Order of A: ' num2str(M.nA)];
   spc = '       ';
   ep.modelEquations = [cr s1 spc s4 cr s2 cr s3 cr];
  elseif ( strcmpi(M.type,'arx') || strcmpi(M.type,'farx') ) 
   cr = sprintf('\n');
   s1 = sprintf('                                   ');
   s2 = sprintf('     A(%s)y(t) = B(%s)u(t) + e(t)',M.op,M.op);
   s3 = sprintf('                                   ');
   s4 = ['Order of A: ' num2str(M.nA)];
   s5 = ['Order of B: ' num2str(M.nB)];
   spc = '       ';
   ep.modelEquations = [cr s1 spc s4 cr s2 cr s3 spc s5 cr];
  end
  

  %----------------------------------------------------------------------
  %  FIR
  %----------------------------------------------------------------------
 case {'fir','nfir'},
  ep.alg='fir';
  for k=1:M.nu,
   if ~strcmpi(M.in(k).type,'linear'),
    ep.startNL = 'startNL';
    ep.startG  = 'startG';    
    ep.alg     = OPT.alg;
   end
  end
  for k=1:M.ny,
   if ~strcmpi(M.out(k).type,'linear'),
    ep.startNL = 'startNL';
    ep.startG  = 'startG';
    ep.alg     = OPT.alg;    
   end
  end
  if M.nu>1,
   g1 = sprintf(' %i               ',M.nu);
   g2 = sprintf('sum  B_i(%s)%s',M.op,inp);
   g3 = sprintf('i=1              ');
  else
   g1 = sprintf('          ',M.op);
   g2 = sprintf('B(%s)%s',M.op,inp);
   g3 = sprintf('          ',M.op);
  end
  h1 =         '         ';
  h2 = sprintf(' + e(t)',M.op);
  h3 =         '         ';
  cr = sprintf('\n');
  s1 = ['            ' g1 h1];
  s2 = [sprintf('      %s = ',out) g2 h2];
  s3 = ['            ' g3 h3];
  s4 = ['Order of B: ' num2str(M.nB(:)')];
  spc = '       ';
  ep.modelEquations = [cr ham cr s1 cr s2 spc s4 cr s3 cr cr wen cr];

  
  %----------------------------------------------------------------------
  %  STATIC
  %----------------------------------------------------------------------
 case {'static'},
  ep.alg=OPT.alg;
  cr = sprintf('\n');
  if M.nu>1,
   g1 = sprintf(' %i               ',M.nu);
   g2 = sprintf('sum %s',M.op,inp);
   g3 = sprintf('i=1              ');
  else
   g1 = sprintf('    ');
   g2 = sprintf('%s',inp);
   g3 = sprintf('    ');
  end
  cr = sprintf('\n');
  s1 = ['            ' g1];
  s2 = [sprintf('     %s = ',out) g2];
  s3 = ['            ' g3];
  ep.modelEquations = [cr ham cr s1 cr s2 cr s3 cr];
  for k=1:M.nu,
   if ~strcmpi(M.in(k).type,'linear'),
    ep.startNL = 'startNL';
   end
  end
  for k=1:M.ny,
   if ~strcmpi(M.out(k).type,'linear'),
    ep.startNL = 'startNL';
   end
  end
  
 %----------------------------------------------------------------------
 %  ARMA, ARMAX, OE, BJ,
 %  (and non-linear versions), NARX, NFIR, NARMA, NARMAX, NOE, NBJ
 %----------------------------------------------------------------------
 case {'arma','armax','oe','bj'},
  
  if ~strcmpi(M.type,'arma'),
   ep.startG='startG';
  end
  ep.startH='startH';
  ep.alg=OPT.alg;
  for k=1:M.nu,
   if ~strcmpi(M.in(k).type,'linear'),
    ep.startNL = 'startNL';
   end
  end
  if ~strcmpi(M.out.type,'linear'),
   ep.startNL = 'startNL';
  end
  
  if strcmpi(M.type,'arma'),
   cr = sprintf('\n');
   s1 = sprintf('                          ');
   s2 = sprintf('     A(%s)y(t) = C(%s)e(t)',M.op,M.op);
   s3 = sprintf('                          ');
   s4 = ['Order of A: ' num2str(M.nA)];
   s5 = ['Order of C: ' num2str(M.nC)];
   spc = '       ';
   ep.modelEquations = [cr s1 spc s4 cr s2 cr s3 spc s5 cr];
  elseif strcmpi(M.type,'armax'),
   if M.nu>1,
    g1 = sprintf(' %i               ',M.nu);
    g2 = sprintf('sum  B_i(%s)%s',M.op,inp);
    g3 = sprintf('i=1              ');
   else
    g1 = sprintf('          ',M.op);
    g2 = sprintf('B(%s)%s',M.op,inp);
    g3 = sprintf('          ',M.op);
   end
   h1 =         '            ';
   h2 = sprintf(' + C(%s)e(t)',M.op);
   h3 =         '            ';
   cr = sprintf('\n');
   s1 = ['            ' g1 h1];
   s2 = [sprintf(' A(%s)%s = ',M.op,out) g2 h2];
   s3 = ['            ' g3 h3];
   s4 = ['Order of B: ' num2str(M.nB(:)') '    Order of C: ' num2str(M.nC)];
   s5 = ['Order of A: ' num2str(M.nA(:)')];
   spc = '       ';
   ep.modelEquations = [cr ham cr s1 spc s4 cr s2 cr s3 spc s5 cr cr wen cr];
  elseif strcmpi(M.type,'oe'),
   if M.nu>1,
    g1 = sprintf(' %i   B_i(%s)       ',M.nu,M.op);
    g2 = sprintf('sum  ------%s',inp);
    g3 = sprintf('i=1  A_i(%s)       ',M.op);
   else
    g1 = sprintf('B(%s)     ',M.op);
    g2 = sprintf('----%s',inp);
    g3 = sprintf('A(%s)     ',M.op);
   end
   h1 = '       ';
   h2 = ' + e(t)';
   h3 = '       ';
   cr = sprintf('\n');
   s1 = ['            ' g1 h1];
   s2 = [sprintf('     %s = ',out) g2 h2];
   s3 = ['            ' g3 h3];
   s4 = ['Order of B: ' num2str(M.nB(:)')];
   s5 = ['Order of A: ' num2str(M.nA(:)')];
   spc = '       ';
   ep.modelEquations = [cr ham cr s1 spc s4 cr s2 cr s3 spc s5 cr cr wen cr];
  elseif strcmpi(M.type,'bj'),
   if M.nu>1,
    g1 = sprintf(' %i   B_i(%s)       ',M.nu,M.op);
    g2 = sprintf('sum  ------%s',inp);
    g3 = sprintf('i=1  A_i(%s)       ',M.op);
   else
    g1 = sprintf('B(%s)     ',M.op);
    g2 = sprintf('----%s',inp);
    g3 = sprintf('A(%s)     ',M.op);
   end
   h1 = sprintf('   C(%s)     ',M.op);
   h2 = ' + ----e(t)';
   h3 = sprintf('   D(%s)     ',M.op);
   cr = sprintf('\n');
   s1 = ['            ' g1 h1];
   s2 = [sprintf('     %s = ',out) g2 h2];
   s3 = ['            ' g3 h3];
   s4 = ['Order of B: ' num2str(M.nB(:)') '    Order of C: ' num2str(M.nC)];
   s5 = ['Order of A: ' num2str(M.nA(:)') '    Order of D: ' num2str(M.nD)];
   spc = '       ';
   ep.modelEquations = [cr ham cr s1 spc s4 cr s2 cr s3 spc s5 cr cr wen cr];
  end
  
 %----------------------------------------------------------------------
 %  State-space and BILINEAR
 %----------------------------------------------------------------------
 case {'ss','bilin','bilinear','lpv'},
  switch OPT.alg,
   %--------------------------------------------------------------
   %  Subsapce
   %--------------------------------------------------------------
   case {'sid','n4sid','cca','subspace'}
    ep.alg='subspace';

    %--------------------------------------------------------------
    %  Guass-Newton or Expectation Maximisation
    %--------------------------------------------------------------
   case {'gn','em'},
    ep.startG='startG';
    ep.startH='startH';
    ep.alg=OPT.alg;
  end

  %Set the model equations
  if M.op=='s',
   s0 = '     u(s)   = u(t)  t <= s < t+d     assume piecewise constant input';
   s1 = '     .                         .   ';
   s2 = '     x(t)   = Ax(t) + Bu(t) + Ke(t)';
   s3 = '     .                         .   ';
   s4 = '     z(t)   = Cx(t)         +  e(t)';
   s5 = '               t+d .';
   s6 = '     y(t+d) =  int z(s)ds            assume integrated sampling ';
   s7 = '                t ';
   cr = sprintf('\n');
   ep.modelEquations = [cr s0 cr cr s1 cr s2 cr cr s3 cr s4 cr cr s5 cr s6 cr s7 cr cr];
  else
   if M.nu>0,
    Bu = ' + Bu(t)';
    if M.estD,
     Du = ' + Du(t)';
    else
     Du = '';
    end
    if M.estF,
     Fukx = ' + F*kron(u(t),x(t))';
    else
     Fukx = '';
    end
    if M.estG,
     Gukx = ' + G*kron(u(t),x(t))';
    else
     Gukx = '';
    end
   else
    Fukx = '';
    Gukx = '';
    Bu   = '';
    Du   = '';
   end
   if M.estK,
    Ke = ' + Ke(t)';
   else
    Ke = '';
   end
   if any(strcmpi(M.type,{'bilin','bilinear'}))
    ep.modelEquations = [sprintf('\n     %sx(t) = Ax(t)',M.op) Bu Fukx Ke '   ' sprintf('nx = %i,  nu = %i,  ny = %i',M.nx,M.nu,M.ny) sprintf('\n      y(t) = Cx(t)') Du Gukx sprintf(' +  e(t)\n',M.op)];
   else
    ep.modelEquations = [sprintf('\n     %sx(t) = Ax(t)',M.op) Bu Ke  '   ' sprintf('nx = %i,  nu = %i,  ny = %i',M.nx,M.nu,M.ny) sprintf('\n      y(t) = Cx(t)') Du sprintf(' +  e(t)\n',M.op)];
   end
  end
   
 %----------------------------------------------------------------------
 %  State-space and BILINEAR
 %----------------------------------------------------------------------
 case {'nlss'},
        
  %Set empty model equations
  ep.modelEquations = []; 
   
  %Make est call the main routine emnlss
  ep.alg     = 'emnlss';
  ep.finishM = 'nlssfinish';
                
  %----------------------------------------------------------------------
  %  Otherwise we do not know the type?
  %----------------------------------------------------------------------
 otherwise,
  error('M.type is unknown!');
end

if ~isfield(ep,'modelEquations'),
 ep.modelEquations = '';
end
