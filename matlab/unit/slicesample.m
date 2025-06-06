%  SLICESAMPLE:  Slice sampler for generating realisations theta_1,
%  theta_2,..... whose distribution converges to an arbitrary density 
%  p(theta|Z) which may be specified by the user.
%
%  Usage is 
%
%  G = slicesample(Z,M,OPT)
%  
%  where:
%
%   Z:          A matlab structure which contains the data which is 
%               used in the conditioning in the density p(theta|Z) that 
%               this routine is seeking to compute. The format is
%               arbitrary, but must be consistent with what the user
%               expects in the user defined function M.pratio
%   M:          A matlab structure which
%

%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%        		      Australia.

%
% Copyright (C) Brett Ninness


function G = slicesample(Z,M,OPT)

% Unspecified parts of OPT -> defaults
if ~exist('OPT') OPT = startOPT([]); else OPT = startOPT(OPT); end;
if ~isfield(OPT,'Mmax')   OPT.Mmax=1e5;               end;
if ~isfield(OPT,'dens')   OPT.dens='gaussian';        end;
if ~isfield(OPT,'mcvar')  OPT.mcvar=mcvar;            end;
if ~isfield(OPT,'burn')   OPT.burn=0.1;               end;

% Set aside memory to store sampler realisation from parametrization of and initialise first column
if isfield(M,'theta')
 theta = M.theta;  thetan = length(theta);
else
 error('Must specify starting value M.theta!');
end;

if isfield(OPT,'var')
 var = OPT.var;  
else
 error('Must specify starting value OPT.var!');
end;

G.TH = zeros(length(theta),OPT.Mmax); G.TH(:,1)=theta(:); thn = length(theta); 
G.varlog = zeros(1,OPT.Mmax);

% Set aside memory to store realisations of noise variance
opt=OPT; G.varlog=zeros(1,OPT.Mmax); G.varlog(1) = OPT.var;

idx =  2;                  % Where we are up to in recording a MC realisation; 
pcom = 0;                  % Percentage complete count initialised to zero;
mark = cputime;            % Used to keep track of elapsed time

% Get target density value at initialisation point
M.thnew = theta; Pstar = feval(M.ptarget,Z,M,OPT);

% Initialise candidate models 
Mleft = M; Mright = M; Mprime = M;
OPTleft = OPT; OPTright = OPT; OPTprime = OPT;

width = 0.01;  % Width for stepping out on each axis.


% Now ready to run the slice sampler
if OPT.dsp
 disp('Running Slice Sampler........')
 disp('');
end;

for k=2:OPT.Mmax  
 
% If requested, give feedback on our status
 if OPT.dsp 
  pcent=10; dcent=100/pcent;
  if ( (mod(k,floor(OPT.Mmax/dcent))==0)|k==2)
   disp(sprintf('Percentage Complete = %d%%, Time since last update = %f s',pcom,cputime-mark)); 
   pcom=pcom+pcent;
   if k>2
    remaining = (dcent - k/floor(OPT.Mmax/dcent))*(cputime-mark);
    hrs  = floor(remaining/3600); remaining = rem(remaining,3600);
    mins = floor(remaining/60); 
    secs = floor(rem(remaining,60));
    disp(sprintf('Predicted completion in %d:%d:%d hrs:mins:secs',hrs,mins,secs))
   end;
   mark=cputime;
  end;
 end; 
 
 % Drop a line a random distance down from Pstar 
 Puprime = Pstar*rand;
 
 Puprime = Pstar + log(rand);
 
 for i=1:thetan   % Sample from each value in parameter vector in turn
  
  Mleft.thnew = theta; Mright.thnew = theta;  Mprime.thnew = theta;
  
  % Generate a random horizonal interval around current value 
  bit = rand;
  Mleft.thnew(i)  = theta(i) - bit*width; 
  Mright.thnew(i) = theta(i) + (1-bit)*width; 
   
  % Now step out until our horizontal interval spans the target density
  while ( feval(M.ptarget,Z,Mleft,OPT) > Puprime)
   Mleft.thnew(i) = Mleft.thnew(i) - width;
  end;

  while ( feval(M.ptarget,Z,Mright,OPT) > Puprime)
   Mright.thnew(i) = Mright.thnew(i) + width;
  end;
  
  stepcount = 0;
  while 1
   stepcount = stepcount+1;
   %fprintf('Iteration %d Step %d     \r',i,stepcount);
   
   % Draw a candidte value uniformly distributed in interval [Mleft,Mright]
   Mprime.thnew(i) = rand()*(Mright.thnew(i) - Mleft.thnew(i)) + Mleft.thnew(i);
   Pstar = feval(M.ptarget,Z,Mprime,OPT);
   
   if (Pstar > Puprime)  % Have we generated a sample under the curve?
    break;  % If so, leave while loop
   else  % If not, reverse the step out - shrink in
    if  (Mprime.thnew(i) > theta(i))
     Mright.thnew(i) = Mprime.thnew(i);
    elseif  (Mprime.thnew(i) < theta(i))
     Mleft.thnew(i) = Mprime.thnew(i);
    else
     error('What?  Cannot find a theta point under target probability');
    end;
   end;
  end;  % While loop that we exit only if we find a Puprime<Pstar
  
  % The random sample found under the target is used to update theta
  theta(i) = Mprime.thnew(i);
  
 end;  % Loop over the i number of terms in theta
 
 % Record this realsation of the whole vector theta
 G.TH(:,k) = theta(:);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % Now do exactly the same for the variance 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 estvar = 1;
 
 if estvar
  % Generate a random horizonal interval around current value 
  bit = rand;
  OPTleft.var  = (sqrt(var) - bit*width)^2; 
  OPTright.var = (sqrt(var) + (1-bit)*width)^2; 
   
  % Now step out until our horizontal interval spans the target density
  while ( feval(M.ptarget,Z,Mprime,OPTleft) > Puprime)
   OPTleft.var = OPTleft.var - width;
  end;
 
  while ( feval(M.ptarget,Z,Mprime,OPTright) > Puprime)
   OPTright.var = OPTright.var + width;
  end;
  
  stepcount = 0;
  while 1
   stepcount = stepcount+1;
   %fprintf('Iteration %d Step %d     \r',i,stepcount);
   
   % Draw a candidate value uniformly distributed in interval [Mleft,Mright]
   OPTprime.var = rand()*(OPTright.var - OPTleft.var) + OPTleft.var;
   Pstar = feval(M.ptarget,Z,Mprime,OPTprime);
   
   if (Pstar > Puprime)  % Have we generated a sample under the curve?
   break;  % If so, leave while loop
   else  % If not, reverse the step out - shrink in
    if  (OPTprime.var > var)
     OPTright.var = OPTprime.var;
    elseif  (OPTprime.var < var)
     OPTleft.var = OPTprime.var;
    else
     error('What?  Cannot find a var point under target probability');
    end;
   end;
  end;  % While loop that we exit only if we find a Puprime<Pstar
  
  % The random sample found under the target is used to update theta
  var = OPTprime.var;  OPT.var = var;
 
  % Record this realsation of the whole vector theta
  G.varlog(k) = var;
 end;  %test on estvar
  
end; % Loop on k up to OPT.Mmax;

G.prop  = 1;
G.mcvar = [];


