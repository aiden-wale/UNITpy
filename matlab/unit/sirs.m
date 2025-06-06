%   SIRS.  This function implements a Sequential Importance Resampling Smoother (SIRS),
%   also known as a particle smoother, for the general model
%
%   x_{t+1} ~ p_1(x_{t+1}|x_t)
%       y_t ~ p_2(y_t|x_t)
%
%   Where p_1 and p_2 are probability density functions. These are
%   defined by a function written by the user, which is passed to this
%   routine via its function handle - see M.model below.
%
%
%   Usage is:  G = sir(Z,M,OPT);
%
%   where
%   
%   Z:           Input-Output data in one of two forms.  The standard form 
%                is for it to be a record with elements Z.y and Z.u, each
%                of which are matrices with number of rows equal to the
%                number of data samples, and number of columns equal (respectively)
%                to the number of outputs and the number of inputs.  On
%                the other hand, Z can be a matrix of the form Z = [y,u] 
%                where it is assumed that y is a column vector of output 
%                measurements and u is a matrix whose columns are the
%                input measurements; in this latter MISO models are 
%                being considered.  
%   M:           Data structure which defines the above model:
%   M.delay:     Number of samples of delay to include (see above model).  
%                In the case of a multi input system, this should be a vector of 
%                delays, one for each input being considered.  Default is all delays 
%                equal to zero.
%   M.model:     Function handle (eg, M.model=@func) for function that
%                specifies the model according to the above distributions
%                p_1 and p_2.  It must be of the form 
%
%                 [q,yp,X2]=func(mode,Z,M,OPT,X,X1)
%
%                where Z,M and OPT are the same structures as passed to
%                SIR (they are passed on the func as a way of passing
%                data Z, model parameters M and optional specs OPT), and
%                X is an array of particles.  Each column represents a
%                particle, and hence it has width equal to the number of
%                particles (OPT.pnum) and height equal to the state
%                dimension.
%  
%
%                The function M.model must conform to the following:
%
%                mode==0:  the function returns the state dimension in 
%                          q (i.e. q = state dimension), all other output
%                          arguments are ignored.
%
%                mode==1:  the returned value q is an array of probabilities:
%
%                          q(i) = p_2(y_t|x_t^i)
%
%                          where x_t^i is specified by the i'th column of X.  
%                          That is, it is the i'th particle at time t.
%
%                          yp is a column vector of expected model output(s) for 
%                          the given particles X, i.e.
%
%                          yp = E{ y_t | X }
%
%                          X2 is ignored in mode 1.                         
%
%                mode==2:  this should have exactly the same functionality as for
%                          mode 1, with the addition that the returned matrix X2 
%                          are sample realisations of particles x_{t+1}^i according 
%                          to the model 
%
%                          x_{t+1} ~ p_1(x_{t+1}|x_t)
%
%                          That is, for each column of X, and with the i'th 
%                          one representing a realisation of the particle x_t^i, 
%                          a new particle x_{t+1}^i is formed by drawing from the
%                          probability distribution 
%
%                          x_{t+1}^i ~ p_1(.|x_t^i).
%
%                mode==3:  the returned value q is an array of
%                          probabilities:
%
%                          q(i) = p_1( x^i_{t+1|N} | x^i_{t|t} )
%
%                          where x^i_{t+1|N} is the i'th column of X1 and
%                          x^i_{t|t} is the i'th colummn of X.
%
%                          yp(:,i) should contain the prediction x^i_{t+1|t}
%
%                mode==4:  the returned variable q is an array of
%                          probabilities:
%
%                          q(i) = p_1( x^i_{t+1|N} | x^i_{t+1|t} )
%
%                          where x^i_{t+1|N} is the i'th column of X1 and
%                          x^i_{t+1|t} is the i'th colummn of X (provided 
%                          by the previous call with mode==3)
%
%                As an example of how to write func, see ssmods.m which
%                implements the LTI model detailed above.
%
%  M.ss.X1:      Initial *predicted* state estimate particles. That is M.ss.X0 =
%                OPT.pnum particles (columns of M.ss.X0) of an initial estimate 
%                of the state at time t=1 given data up to t=0.  
%                Default is M.ss.X0 = zeros(nx,OPT.pnum); where nx is the state
%                dimension.
%
%   OPT:         Data structure which defines options for the estimation
%                algorithm as follows:
%    OPT.pnum    Number of particles to use.  Default is OPT.pnum=100.
%
%    OPT.allP    If nonzero, then the output structure G contains an element
%                G.Pt which is a cell array containing the history of
%                the filtered state sequence covariance matrix.
%
%   G:           Data structure which specifies the Kalman Filter Output as follows:
%                
% G.ss.Xp:       Time evolution of one-step ahead *predicted* state sequence estimates.
% G.ss.Xf:       Time evolution of filtered state sequence estimates.
% G.ss.Xs:       Time evolution of smoothed state sequence estimates.
%   G.yp:        One step ahead measurement prediction y_{t|t-1}
%   G.yf:        Filtered measurement estimate y_{t|t}.
%   G.ys:        Smoothed measurement estimate y_{t|N}.
%
%                             

% Copyright (C) Brett Ninness

function G = sirs(Z,M,OPT);

% Extract sizes of input and output from data matrix
Z = startZ(Z);
[y,u,ny,nu,Ny] = Z2data(Z);

% Model defined by passing a function handle?  If not, assume LTI ss description in M
if ~isfield(M,'model') error('Need M.model to be specified!');           end;

% Set some generic defaults so we play nicely with other toolbox children
if ~isfield(M,'op')     M.op = 'q';                                     end;
if ~isfield(M,'T')      M.T=1;                                          end;   
if ~isfield(M,'delay')  M.delay=zeros(nu,1);                            end; 

% Check what algorithm options are not specified explicitly and set to defaults
if ~exist('OPT')
  OPT.pnum = 100;          % Number of particles
  OPT.allP = 0;            % Don't save all cov matrices as default
else 
  if ~isfield(OPT,'pnum') OPT.pnum = 100;      end;   
  if ~isfield(OPT,'allP') OPT.allP = 0;        end;   
end;

% Copy model specification into output
G=M;

% Include delays specified in model structure on inputs
for r=1:nu u(:,r) = [zeros(M.delay(r),1);u(1:Ny-M.delay(r),r)]; end

%Get state dimension
nx = feval(M.model,0,Z,M);

%Set some local variables
N=Ny; MM=OPT.pnum;

%Get initial state particles
X = G.ss.X1;

%Initialise index
idx=uint32(zeros(OPT.pnum,1));
syst_res=0:OPT.pnum-1;

% Make some room for things (good test for memory requirements)
G.ss.Xp=zeros(nx,N+1);        % Predicted state
G.ss.Xf=zeros(nx,N+1);        % Filtered state 
G.ss.Xs=zeros(nx,N+1);        % Smoothed state
G.ss.PXf=zeros(nx,MM,N+1);     % Filtered particles
G.ss.PXs=zeros(nx,MM,N+1);     % Smoothed particles
q=zeros(MM,1);                 % Weights

%---------------------------------------------
%  RUN FILTER
%---------------------------------------------
% Extract initial predicted state
G.ss.Xp(:,1)=mean(X,2);

G.LL=0;

%Determine if matlab or C based resampling is to be used
if exist('resampling')==3,
 c_or_mat='c';
else
 c_or_mat='m';
end

% Iterate through data updating particles
for t=1:Ny
    % Compute Measurement update weights q_t^i = p(y_t|x^i_t|t-i)  
    z.y=y(t,:); z.u=u(t,:); 
    [q,G.yp(:,t),PY]=feval(M.model,1,z,M,OPT,X);
  
 %Compute LL
 G.LL=G.LL-0.5*(log(det(PY))+(Z.y(t,:)'-G.yp(:,t))'*(PY\(Z.y(t,:)'-G.yp(:,t))));
 
    % Systematic resampling
 % Resample particles on last iteration according to density implied by q
 rspl=1;
 switch rspl,
  case 0,
   % Simple resampling
   ut=cumprod(rand(1,OPT.pnum).^(1./(OPT.pnum:-1:1))); ut=fliplr(ut);
  case 1,
   % Systematic resampling
   ut=(syst_res+rand)/OPT.pnum;
  case 2,
   % Brett's resampling
   ut=rand(1,OPT.pnum);
 end
 cp=cumsum(q);
 switch c_or_mat,
  case 'c',
   resampling(cp,ut,idx);
  case 'm'
   i=1;for k=1:OPT.pnum,while cp(i)<ut(k),i=i+1;end;idx(k)=i;end
 end
    X=X(:,idx); % Resample from particles at time t-1 to get particles at time t
 
    % Extract filtered state
    G.ss.Xf(:,t)=mean(X,2); G.ss.PXf(:,:,t)=X;
 
    % Now do time update of particles x^i_t+1|t \sim p(x_t+1|x_t^i) 
    [q,G.yf(:,t),X]=feval(M.model,2,z,M,OPT,X);

    % Extract predicted state
    G.ss.Xp(:,t+1)=mean(X,2);
end;
G.ss.Xf(:,t+1)=mean(X,2); G.ss.PXf(:,:,t+1)=X;


%---------------------------------------------
%  RUN SMOOTHER
%---------------------------------------------
% Set initial conditions
G.ss.Xs(:,N)=G.ss.Xf(:,N);
G.ss.Xs(:,N+1)=G.ss.Xf(:,N+1); 
G.ss.PXs(:,:,N)=G.ss.PXf(:,:,N);
G.ss.PXs(:,:,N+1)=G.ss.PXf(:,:,N+1); 

% Run backwards filter
for t=N-1:-1:1,
    % Get Particles for x_{t+1|N}, x_{t|t}
    xN=squeeze(G.ss.PXs(:,:,t+1)); 
    xt=squeeze(G.ss.PXf(:,:,t)); 
    
    % p_1(x^i_{t+1|N} | x^i_{t+1|t})
    z.u=u(t,:);
    [p,xp] = feval(M.model,3,z,M,OPT,xt,xN);
    
    for i=1:MM,
        % return vector of p_1(x_{t+1|N}^i | xp^i) for all i
        num=(1/MM)*sum(feval(M.model,4,z,M,OPT,xp,repmat(xN(:,i),1,MM)));
        q(i)=p(i)/num;
    end;

    % Systematic resampling
 % Resample particles on last iteration according to density implied by q
 rspl=1;
 switch rspl,
  case 0,
   % Simple resampling
   ut=cumprod(rand(1,OPT.pnum).^(1./(OPT.pnum:-1:1))); ut=fliplr(ut);
  case 1,
   % Systematic resampling
   ut=(syst_res+rand)/OPT.pnum;
  case 2,
   % Brett's resampling
   ut=rand(1,OPT.pnum);
 end
 cp=cumsum(q);
 switch c_or_mat,
  case 'c',
   resampling(cp,ut,idx);
  case 'm'
   i=1;for k=1:OPT.pnum,while cp(i)<ut(k),i=i+1;end;idx(k)=i;end
 end
 G.ss.PXs(:,:,t)=G.ss.PXf(:,idx,t); 
 G.ss.Xs(:,t)=mean(G.ss.PXs(:,:,t),2);
end

% Work out smoothed measurement estimate

G.ys = M.ss.C*G.ss.Xs(:,2:end) + M.ss.D*u';



