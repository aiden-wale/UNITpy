%  METROPOLIS:  Metropolis algorithm implementation for generating a
%  Markov Chain theta_1, theta_2,..... whose distribution converges to 
%  an arbitrary density p(theta|Z) which may be specified by the user.
%
%  Usage is 
%
%  G = metropolis(Z,M,OPT)
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
%        		              Australia.

%
% Copyright (C) Brett Ninness

function G = metropolis(Z,M,OPT)

mcvar = 1e-6;        % default variance of random walking driving MC.
Temp = 1000;         % Initial temperature for annealing
dfac = 0.92;         % Factor to decrease temperature by
alfwin=100;          % Width of sliding window used to monitor average acceptance rate

% Unspecified parts of OPT -> defaults
if ~exist('OPT') OPT = startOPT([]); else OPT = startOPT(OPT); end;
if ~isfield(OPT,'Mmax')   OPT.Mmax=1e5;               end;
if ~isfield(OPT,'dens')   OPT.dens='gaussian';        end;
if ~isfield(OPT,'mcvar')  OPT.mcvar=mcvar;            end;
if ~isfield(OPT,'burn')   OPT.burn=0.1;               end;
if ~isfield(OPT,'rej')    OPT.rej=0;               end;
OPT.cold=[];  % Will be used as a way of telling M.pratio not to bother computing old cost

% Set aside memory to store chain realisation and initialise first column
if isfield(M,'theta')
 theta = M.theta;
else
 error('Must specify starting value M.theta!');
end;
G.TH = zeros(length(theta),OPT.Mmax); G.TH(:,1)=theta(:); thn = length(theta); 
G.mcvar = zeros(1,OPT.Mmax); % Record of (possibly) dynamic changing proposal variability

% To estimate noise variance, keep metropolis realisations in opt.var
opt=OPT; G.varlog=zeros(1,OPT.Mmax); G.varlog(1) = OPT.var;

idx =  2;                  % Where we are up to in recording a MC realisation; 
pcom = 0;                  % Percentage complete count initialised to zero;
qvar = OPT.mcvar;          % Starting proposal density variance;
accep = 0;                 % Initialize count of number of acceptances.
wincount = 0;              % Acceptance count just over analysis window 
P2=rchol(M.P);             % Scaling on parameter updates
mark = cputime;            % Used to keep track of elapsed time

% Now ready to run the chain
if OPT.dsp
 disp('Running Metropolis Sampler........')
 disp('');
end;
for k=2:OPT.Mmax  
 % If requested, give feedback on our status
 if OPT.dsp 
  if ( (mod(k,floor(OPT.Mmax/20))==0)|k==2)
   disp(sprintf('Percentage Complete = %d%%, Time since last update = %f s',pcom,cputime-mark)); pcom=pcom+5;
   if k>2
    remaining = (20 - k/floor(OPT.Mmax/20))*(cputime-mark);
    hrs  = floor(remaining/3600); remaining = rem(remaining,3600);
    mins = floor(remaining/60); 
    secs = floor(rem(remaining,60));
    disp(sprintf('Predicted completion in %d:%d:%d hrs:mins:secs',hrs,mins,secs))
   end;
   mark=cputime;
  end;
 end; 

 % Every time a block of width alfwin passes, re-examine proposal variance
 if ~mod(k,alfwin)  
  prop = wincount/alfwin;
  if prop < 0.25
   qvar = qvar/1.2;
  elseif prop > 0.20
   qvar = qvar*1.2;
  end;
  % Gradually increase the width alfwin used for monitoring acceptance rate 
  %if k>floor(0.1*OPT.Mmax)
  if k>10000 
   alfwin=2*alfwin;
  end;
  wincount=0;  % Start counting over again on new window;
 end; 
 
 G.mcvar(k) = qvar; % Keep record of the variances for possible later examination
 
 % Set up mean and variance for proposal
 M.thold = G.TH(:,k-1);  P=(diag(abs(M.thold)));
  
 % Draw new proposal for parametrization of dynamics 
 %M.thnew = M.thold + qvar*randn(size(M.thold));    
 M.thnew = M.thold + qvar*P2'*randn(size(M.thold));   
 
 % Draw new proposal for measurement noise variance 
 varold = opt.var;
 opt.var = (sqrt(opt.var)+sqrt(0.0001)*randn)^2;   
 opt.varold = varold;
 
 [prat,cold]=feval(M.pratio,Z,M,opt); 
 
 % Now keep the new thnew with the associated Metropolis acceptance probability
 if (rand <= prat) 
  G.TH(:,k) = M.thnew;   % Succesful in testing Uniform against alpha
  accep=accep+1;         % Global acceptance count
  wincount=wincount+1;   % Local window acceptance count
  OPT.cold=[];           % Tell M.pratio it needs to recompute "old" cost next time
  G.varlog(k) = opt.var; % Record that we accept this new noise variance  
 else 
  G.TH(:,k) = M.thold;   % The converse, we made a rejection
  opt.var = varold;      % New noise variance was not accepted
  G.varlog(k) = varold;
  OPT.cold = cold;       % No acceptance means denominator of alpha remains unchanged - tell M.pratio no need to recompute it
 end;

end; % Loop on k up to OPT.Mmax;

G.prop = accep/OPT.Mmax;  % Return proportion of accepted proposals


