%  POSTDIST: Compute the posterior distribution marginals p(theta_k|y) of parameters
%  in a dynamic model structure given observed data.  Importantly, it is
%  also possible to compute the posterior distribution p(g(theta)|y) of
%  an arbitrary function g(theta) of the parameters.  The underlying
%  method is Monte-Carlo based via the Metropolis-Hastings algorithm, and the 
%  supported model structure is
%
%          B(p)                   C(p) 
%  y_t  =  ---- u_{t-delay}   +   ---- e_t
%          A(q)                   D(p)
%
%  Usage is 
%
%  G = postdist(Z,M,OPT)
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
%               description.  Typically, the model structure used here
%               would come from a preceding call to est.m in order to
%               estimate a model - see the usage example below for an
%               illustration of this.
%
%  OPT:         Data structure which defines options for the underlying
%               Metropolis-Hastings algorithm.
%  OPT.Mmax:    Mmax - the number of iterations of the Markov Chain
%               implemented by the Metropolis Hastings method from which
%               sample histograms form the posterior densities computed
%               by this postdist.m.  Default is OPT.Mmax=1e5; 
%  OPT.dens:    Specification of the density for the measurement noise 
%               innovations.  Possibilitie are
%               OPT.dens = 'gaussian' - self explanatory, and the default;
%               OPT.dens = 'uniform' - Uniform and zero mean;
%  OPT.mcvar:   Initial variance for proposal density q in Metropolis
%               Hastings Algorithm.  This is subsequently refined by the
%               algorithm to aim at a 55%-65% acceptance rate.
%  OPT.plot:    A flag, that if set will cause the computed posterior 
%               marginals to be plotted, together with corresponding 
%               Normal approximation resulting from usual asymptotic analysis.
%               Default is OPT.plot=0;  Note that plots can alternatively
%               be produced with plotdist(p) call subsequent to running 
%               p = postdist(z,m);
%  OPT.burn:    Proportion of samples to be considered part of the "burn
%               in" period and hence to be thrown away.  Default is 0.2
%  G:           Data structure which specifies the estimated posterior densities
%               as follows.  Note that G will also contain all the
%               elements passed in by the model structure M as well.
%  G.pa(i).p    Posterior marginal density for the i'th element of the denominator
%               estimate A(q).  Note thtat G.pb(i).p, G.pc(i).p and P.pd(i).p are also 
%               provided (if applicable) for other parameters in the model
%  G.pa(i).x    Corresponding x axis for the above.  That is
%               plot(G.pa(1).x.G.pa(1).p) provides a plot of the posterior marginal 
%               p(a_1|Y).
%  G.pex(i).x   Same as above, but "exact" results computed by numerical 
%  G.pex(i).p   integration of joint density to give marginals - provided
%               OPT.int is non-zero.
%
%
%   Usage Example:     
%
%               z.y=y; z.u=u; m.A=1;  % Specify data and model structure
%               g=est(z,m);           % Estimate a 4th order model
%               p=postdist(z,g);      % Compute posterior dist of parameters
%                                     % and display them.  
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%        		              Australia.

%
% Copyright (C) Brett Ninness


function G = postdist(Z,M,OPT)

mmax  = 1e4;         % Default number of runs of chain
dens  = 'gaussian';  % Default density assumed for measurement noise
mcvar = 1e-4;        % default variance of random walking driving MC.
Temp = 1000;         % Initial temperature for annealing
dfac = 0.92;         % Factor to decrease temperature by

% Extract sizes of input and output from data matrix
[y,u,ny,nu,Ny] = Z2data(Z);

% Unspecified parts of OPT -> defaults
if ~exist('OPT') OPT = startOPT([]); else OPT = startOPT(OPT); end;
if (OPT.n>=Ny) error('Cannot have OPT.n larger than number of data samples!'); end;
if ~isfield(OPT,'Mmax')    OPT.Mmax=mmax;              end;
if ~isfield(OPT,'dens')    OPT.dens='gaussian';        end;
if ~isfield(OPT,'mcvar')   OPT.mcvar=mcvar;            end;
if ~isfield(OPT,'plot')    OPT.plot=0;                 end;
if ~isfield(OPT,'int')     OPT.int=0;                  end;
if ~isfield(OPT,'burn')    OPT.burn=0.2;               end;
if ~isfield(OPT,'sampler') OPT.sampler='metropolis';   end;

% Unspecified parts of M -> defaults
if ~exist('M') error('Need to specify initial model structure M!'); 
else M = startM(Z,M);  end;

% Check to see if only integer orders were specified as initial guesses
% for dynamics: if so get initial estimate by fitting ARX model structure.
M = startG(Z,M,OPT);

% Check to see of only integer orders where specified as initial guesses
% for noise model: if so get initial estimate via Hannan-Rissanen method.
M = startH(Z,M,OPT);

% Take parameters in initial model structure and stack them into a parameter vector. 
theta = m2theta(M);   M.theta=theta; 

% Unspecified parts of regularisation model -> defaults
if isfield(OPT,'M') OPT.M = startM(Z,OPT.M); else OPT.M = theta2m(theta*0,M); end; 
OPT.step = 1; OO = OPT; OO.dsp = 0;

% If variance of noise not specified, estimate from residuals;
if ~isfield(OPT,'var') g = est(Z,M,OO); OPT.var=g.var; end;

% Call Metropolis Algorithm to compute posterior

t=cputime;       % Start stopwatch to time how long metropolis takes to run  
if strcmpi(OPT.sampler,'slice')  % Running a slice sampler
 M.ptarget = @ptarget;  % Specify subroutine defining target density
 gm = slicesample(Z,M,OPT);
else                    % Otherwise run a Metropolis sampler
 M.pratio = @pratio;    % Specify subroutine defining target prob ratio. 
 gm = metropolis(Z,M,OPT); 
end;
tm = cputime-t;  % Stop stopwatch to time how long metropolis takes to run  
%gm = rejection(Z,M,OPT);

% Throw away samples from burn in period
[dummy,len] = size(gm.TH); len = floor(len*OPT.burn);  % Compute number of samples to throw away
gm.TH = gm.TH(:,len+1:end);

% Ask theta2m to add marginal posterior density estimates to G
[g,sd] = theta2m(gm.TH,M);

% If requested, plot the results.
if OPT.plot showdist(g); end;  

% Compute parameter point estimate as posterior mean
th = mean(gm.TH');  th=th(:); % Posterior mean estimate in parameter space
G = g;                        % Fill in all incoming info about model structure
G = theta2m(th,G,0);          % Now hand back transfer function version of posterior mean
G = m2f(G);                   % Finally, hand back associated frequency response
G.disp.legend = 'Posterior Mean Estimate';

% Pass back other relevant information 
 
G.TH = gm.TH;         % Markov Chain realisation from metropolis step.
G.prop = gm.prop;     % Proportion of actual acceptances
G.mcvar = gm.mcvar;   % Time history of adaptive proposal variance
G.cputime.tm = tm;    % Record time taken to run Metropolis algorithm

G.varlog = gm.varlog;