%   SIR.  This function implements a Sequential Importance Resampler (SIR),
%   also known as a particle filter, for the general model
%
%   x_{t+1} ~ p_1(x_{t+1}|x_t)
%       y_t ~ p_2(y_t|x_t)
%
%   Where p_1 and p_2 are probability density functions. These are
%   defined by a function written by the user, which is passed to this
%   routine via its function handle - see M.model below.
%
%   Alternatively, if M.model is not specified, then this routine uses
%   the linear time invariant model
%
%    x_{t+1} = Ax_t + Bu_(t-delay) + Lw_t
%    y_t = Cx_t + Du_(t-delay) + v_t
%
%   Where w_t and v_t are zero mean vector white Gaussian processes with covariances
%
%   E(w_tw_t^T) = Q, E(v_tv_t^T) = R, E(w_tv_t^T) = S
%
%   and y_t, u_t are observed output and input processes.
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
%                 [q,X1,yp]=func(Z,M,OPT,X)
%
%                where Z,M and OPT are the same structures as passed to
%                SIR (they are passed on to func as a way of passing
%                data Z, model parameters M and optional specs OPT), and
%                X is an array of particles.  Each column represents a
%                particle, and hence it has width equal to the number of
%                particles (OPT.pnum) and height equal to the state
%                dimension.
%
%                The returned value q is an array of probabilities:
%
%                q(i) = p_2(y_t|x_t^i)
%
%                where x_t^i is specified by the i'th column of X.  That
%                is, it is the i'th particle at time t.
%
%                The returned matrix X1 contains sample realisations of
%                particles x_{t+1}^i according to the model
%
%                x_{t+1} ~ p_1(x_{t+1}|x_t)
%
%                That is, for each column of X, and with the i'th one
%                representing a realisation of the particle x_t^i, a new
%                particle x_{t+1}^i is formed by drawing from the
%                probability distribution
%
%                x_{t+1}^i ~ p_1(|x_t^i).
%
%                Finally, yp is a vector of model outputs, with the i'th
%                column being the value associated with the i'th column
%                of X, and hence the i'th particle realisation x_t^i.
%
%                THIS IS IMPORTANT.  In order to inform SIR what the
%                underlying state dimension is, the function specifying
%                the model must check to see if it is called as func(Z,M)
%                - that is, with only two parameters - and in this
%                special case it must return an integer which represents
%                the state dimension.
%
%                As an example of how to write func, see ssmod.m which
%                implements the LTI model detailed above.
%
%                The following elements M.ss are used only if M.model is
%                not specified (or it is set as M.model=@ssmod).
%
%M.ss.A,M.ss.B:  [A,B,C,D] matrices/vectors defining state space description.
%M.ss.C,M.ss.D:  of underlying dynamic system as shown above.
%  M.ss.L:       Matrix defining how state noise enters.  Default is M.ss.L=I.
%M.ss.Q,M.ss.R:  Covariance matrices for state and measurement noise.
%  M.ss.S:       as described above.  Defaults are M.ss.R = I, M.ss.Q=0.01*I, M.ss.S=0;
%  M.ss.X0:      Initial *predicted* state estimate.  That is M.ss.X0 =
%                an initial estimate of state at time t=1 given data up
%                to t=0.  Default is M.ss.X0 = 0;
%  M.ss.P0:      Covariance in initial *predicted* state estimate (see above).
%                Default is M.ss.P0 = 100*I.
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
% G.ss.P:        Final covariance matrix of predicted state estimates.
% G.ss.Pf:       Final covariance matrix of filtered state estimates.
% G.ss.Pt:       If OPT.allP is nonzero this is a cell array containing the history of
%                the filtered state sequence covariance matrix.  That is
%                G.Pt{k} is the covariance P_{k|k}.
%   G.yp:        One step ahead measurement prediction y_{t|t-1}
%   G.yf:        Filtered measurement estimate y_{t|t}.
%   G.pe:        One step ahead output prediction error.
%   G.LL:        Log-Likelihood cost of given data under the model specified
%                by M.
%   G.mse:       Mean square cost of data under the model specified by M.
%
%
% Written by Brett Ninness, School of EE & CS
%                           University of Newcastle
%                           Australia.

% Copyright (C) Brett Ninness.

function G = sir(Z,M,OPT);

% Extract sizes of input and output from data matrix
Z = startZ(Z);
[y,u,ny,nu,Ny] = Z2data(Z);

% Model defined by passing a function handle?  If not, assume LTI ss description in M
if isfield(M,'model') ss=0; else ss=1; end;

% Set some generic defaults so we play nicely with other toolbox children
if ~isfield(M,'op')     M.op = 'q';                                    end;
if ~isfield(M,'T')      M.T=1;                                         end;
if ~isfield(M,'delay')  M.delay=zeros(nu,1);                           end;

% If default ss model, set unspecified bits to defaults
if ~ss
 if ~isfield(M.ss,'A')   error('Need to specify M.ss.A!');              end;
 if ~isfield(M.ss,'B')   error('Need to specify M.ss.B!');              end;
 if ~isfield(M.ss,'C')   error('Need to specify M.ss.C!');              end;
 if ~isfield(M.ss,'D')   M.ss.D = zeros(ny,ny);                         end;
 if ~isfield(M.ss,'R')   M.ss.R = eye(ny,ny);                           end;
 [m,n] = size(M.ss.A);   if ~isfield(M.ss,'Q') M.ss.Q = 0.001*eye(n,n); end;
 [mm,nn] = size(M.ss.Q); if ~isfield(M.ss,'L') M.ss.L = eye(n,nn);      end;
 if ~isfield(M.ss,'S')   M.ss.S = zeros(n,ny);                          end;
 if ~isfield(M.ss,'X0')  M.ss.X0 = zeros(n,1);                          end;
 if ~isfield(M.ss,'P0')  M.ss.P0 = 100*eye(n,n);                        end;
 if ~isfield(M,'type')   M.type='ss';                                   end;
 % Fix up idiotic matlab "feature" that x+[] = [] and not x!
 if isempty(M.ss.D) M.ss.D = zeros(ny,max(nu,1));       end;
 if isempty(M.ss.B) M.ss.B = zeros(size(M.ss.A,1),1);   end;
 if isempty(M.ss.R) M.ss.R = eye(ny);                   end;
end;

% Check what algorithm options are not specified explicitly and set to defaults
if ~exist('OPT')
 OPT.pnum = 100;          % Number of particles
 OPT.allP = 0;            % Don't save all cov matrices as default
 OPT.hist = 0;            % No saving of histograms as default
else
 if ~isfield(OPT,'pnum') OPT.pnum = 100;      end;
 if ~isfield(OPT,'allP') OPT.allP = 0;        end;
 if ~isfield(OPT,'hist') OPT.hist = 0;        end;
end;

% Copy model specification into output
G=M;

% Include delays specified in model structure on inputs
for r=1:nu u(:,r) = [zeros(M.delay(r),1);u(1:Ny-M.delay(r),r)]; end

if ss % Case of default ss model being used (M.model was not defined)
 nx = size(M.ss.A,1);     % State Dimension
else % Ask M.model what the state dimension is
 nx = feval(M.model,Z,M);
end;

%X=zeros(nx,OPT.pnum);    % Initialise particles (all at the origin)
X=repmat(M.ss.X0,1,OPT.pnum)+chol(M.ss.P0)*randn(nx,OPT.pnum);
chol_Q=sqrtm(M.ss.Q);  %Form Cholesky factor outside main loop
idx=uint32(zeros(OPT.pnum,1)); %initialise index
syst_res=0:OPT.pnum-1;

G.LL=0;

%Determine if matlab or C based resampling is to be used
if exist('resampling')==3,
 c_or_mat='c';
else
 c_or_mat='m';
end

% Now iterate through data updating particles
for t=1:Ny
 % Compute Measurement update weights q_t^i = p(y_t|x^i_t|t-i)
 if ss  % Case of default ss model being used (M.model was not defined)
  yp = M.ss.C*X(:,:)+M.ss.D*repmat(u(t,:)',1,OPT.pnum);
  G.yp(:,t)=mean(yp,2);
  pe=repmat(y(t,:)',1,OPT.pnum)-yp;
  qpe = pe.*(inv(M.ss.R)*pe);
  if size(qpe,1)>1 qpe=sum(qpe); end;
  q=exp(-qpe/2); q=q/sum(q);
 else % Use model defined by function handle M.model
  z.y=y(t,:); z.u=u(t,:);
  [q,dummy,G.yp(:,t),PY]=feval(M.model,z,M,OPT,X);
  %Compute LL
  G.LL=G.LL-0.5*(log(det(PY))+(Z.y(t,:)'-G.yp(:,t))'*(PY\(Z.y(t,:)'-G.yp(:,t))));
 end;

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
 X = X(:,idx); % Resample from particles at time t-1 to get particles at time t

 % If user asks for histogram, give it to them
 if OPT.hist G.ss.Xhist(:,:,t) = X; end;

 % Extract filtered state and it's covariance
 G.ss.Xf(:,t)=mean(X,2); G.ss.PXf(:,:,t)=X;
 if t==Ny G.ss.Pf = cov(X',1);        end;
 if OPT.allP G.ss.Pt{t} = cov(X',1);  end;

 % Now do time update of particles x^i_t+1|t \sim p(x_t+1|x_t^i)
 if ss % Case of default ss model being used (M.model was not defined)
  % First extract filtered output estimate
  yf = M.ss.C*X(:,:)+M.ss.D*repmat(u(t,:)',1,OPT.pnum);
  G.yf(:,t)=mean(yf,2);
  X = M.ss.A*X+M.ss.B*repmat(u(t,:)',1,OPT.pnum)+chol_Q*randn(nx,OPT.pnum);
 else % Use model defined by function handle M.model
  [q,X,G.yf(:,t)]=feval(M.model,z,M,OPT,X);
 end;

 % Extract predicted state and its covariance
 G.ss.Xp(:,t+1)=mean(X,2);
 if t==Ny G.ss.P = cov(X',1); end;
end;

G.ss.PXf(:,:,t+1)=X;

%toc

% Fill in remaining bits of output
G.pe = y'-G.yp;
G.mse = (G.pe*G.pe')/(Ny*ny);

