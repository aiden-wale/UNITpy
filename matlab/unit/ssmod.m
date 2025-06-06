%   Function to compute the probability updates
%
%   x_{t+1} ~ p_1(x_{t+1}|x_t)
%       y_t ~ p_2(y_t|x_t)
%
%   for the linear time invariant model
%
%    x_{t+1} = Ax_t + Bu_(t-delay) + Lw_t
%    y_t = Cx_t + Du_(t-delay) + v_t
%
%   Where w_t and v_t are zero mean vector white Gaussian processes with covariances
%
%   E(w_tw_t^T) = Q, E(v_tv_t^T) = R, E(w_tv_t^T) = S
%
%   and y_t, u_t are observed output and input processes.  This function
%   is intended for possible use with the SIR.m sequential importance
%   resampler routine.
%
%   Usage is:  [q,X1,yp]=ssmod(Z,M,OPT,X)
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
%   OPT:         Data structure which defines options for the estimation
%                algorithm as follows:
%    OPT.pnum    Number of particles to use.  Default is OPT.pnum=100.
%
%  X:            Array of particles. Each column represents a
%                particle, and hence it has width equal to the number of
%                particles (OPT.pnum) and height equal to the state
%                dimension.
%  
%  q:            Aarray of probabilities:
%
%                q(i) = p_2(y_t|x_t^i)
%
%                where x_t^i is specified by the i'th column of X.  That
%                is, it is the i'th particle at time t.
%
%  X1:           Matrix containing sample realisations of
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
%  yp:           Vector of model outputs, with the i'th
%                column being the value associated with the i'th column
%                of X, and hence the i'th particle realisation x_t^i.
%
%  IMPORTANT: In order to inform SIR what the underlying state dimension is,
%                this function checks to see if it is called as ssmod(Z,M) -
%                that is, with only two parameters - and in this special
%                case it must return an integer which represents the state
%                dimension.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness, Adrian Wills

function [q,X,yp,PY]=ssmod(Z,M,OPT,X)
 
nx = size(M.ss.A,1);  % Determine state dimension

chol_Q=sqrtm(M.ss.Q); % Form Cholesky factor outside main loop

if nargin<3 % Nargin<3 => user just wants state dimension returned
 q=nx;
else
 % Compute Measurement update weights q_t^i = p(y_t|x^i_t|t-i) 
 yp = M.ss.C*X(:,:)+M.ss.D*repmat(Z.u(:),1,OPT.pnum);
 pe=repmat(Z.y(:),1,OPT.pnum)-yp; 
 
 PY=pe*pe'/OPT.pnum;
 
 qpe = pe.*(M.ss.R\pe);
 if size(qpe,1)>1 qpe=sum(qpe); end;
 q=exp(-qpe/2); q=q/sum(q);        
 
 % Pass back predicted y based on passed in (predicted) state
 yp=mean(yp,2);
 
 % Now do time update of particles x^i_t+1|t \sim p(x_t+1|x_t^i)
 X = M.ss.A*X+M.ss.B*repmat(Z.u(:),1,OPT.pnum)+chol_Q*randn(nx,OPT.pnum); 
end;    
    
