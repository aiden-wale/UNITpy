%   RECUR This routine runs a recursive estimation algorithm on a linear
%   predictor type model.  That is, for the model structure
% 
%         y_t = phi_t^T*theta
%
%   this routine runs a recursive (in time) altorithm for the estimation
%   of theta.  The routine here is not meant to be called directly by the
%   user.  Rather, it is meant as a subroutine of other functions (onid,
%   barx) for situations in which a user requests a recursive solution.
%
%   G = recur(y,PHI,OPT); 
%
%   where
%
%   y    =   is a column vector of (scalar) output measurements.
%   PHI  =   matrix with columns containing the time evolution of the
%            regressor elements in phi_t.  That is, PHI(t,:) = phi_t^T.
%   OPT  =   Data structure which defines options for the estimation
%            algorithm as follows:
%   OPT.n  = number of starting data points to ignore.  Default = 0.
%   OPT.alg.type - determines the algorithm type used.  It may set 
%            to any of:
%            'block' - the solution for the least squares 
%                      estimate is found in block form (default).
%            'rls'   - the solution is found recursively 
%                      via the recursive least squares algorithm. 
%            'ct'    - the solution is found recursively via recursive 
%                      least squares with a contant trace forced
%                      on the covariance matrix.
%            'lms'   - the solution is found recursively via 
%                      least-mean-square (gradient descent) algorithm.
%            'kf'    - the solution is found recursively via 
%                      a Kalman Filtering algorithm.
%   OPT.alg.P      - Initialisation of covariance matrix, default is P=10*I
%   OPT.alg.mu     - LMS gain, default is mu=0.001;
%   OPT.alg.lambda - RLS alg. `forgetting factor'. Default is lambda = 1;
%   OPT.alg.R      - Measurement noise variance used by kf alg. Default = 1;
%   OPT.alg.Q      - Parameter variance used by kf alg. Default = 0.1*I;
%   OPT.alg.th     - Initial parameter vector value.  Default = 0.
%
%  G          = Data structure specifying estimated model as follows:
%    G.th     = Final parameter vectorestimate.
%    G.th_hist= History of how parameter space estimates evolved.
%    G.P      = Final covariance Matrix for Estimated Parameters.
%    G.pe     = Time history of how prediction error evolved.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%        		      Australia.

% Copyright (C) Brett Ninness.
 
function G = recur(y,PHI,OPT); 
 
% Get data length and dimentions of regressor matrix 
y = y(:); [Ny,nin] = size(y); 
[Npsi,d] = size(PHI); if d>Npsi PHI=PHI'; [Npsi,d] = size(PHI); end;

% Unspecified parts of OPT -> defaults
if ~exist('OPT') OPT = startOPT([]); else OPT = startOPT(OPT); end;
if (OPT.n>=Ny) error('Cannot have OPT.n larger than height of Z!'); end;
if ~isfield(OPT,'alg')        OPT.alg.type = 'rls';    end;
if ~isfield(OPT.alg,'P')      OPT.alg.P = 10*eye(d);   end;
if ~isfield(OPT.alg,'mu')     OPT.alg.mu = 0.001;      end;  
if ~isfield(OPT.alg,'th')     OPT.alg.th = zeros(d,1); end;  
if ~isfield(OPT.alg,'lambda') OPT.alg.lambda = 1;      end;  
if ~isfield(OPT.alg,'Q')      OPT.alg.Q = eye(d);      end;  
if ~isfield(OPT.alg,'R')      OPT.alg.R = 1;           end;  

% Initialise algorithm with specs given in OPT structure
P = OPT.alg.P; th = OPT.alg.th;

% Reserve some memory for the saved outputs
G.th_hist = zeros(Ny,d); G.pe = zeros(1,Ny);

% Now run through the data finding the recursive solution.
for t=OPT.n+1:Ny
 phi = PHI(t,:)';                    % Extract regressor for this time step
 e = y(t) - phi'*th;                 % Prediction error at this time step

 if ~strcmp(lower(OPT.alg.type),'lms')
  Pphi  = P*phi;  denom = phi'*Pphi; % Common term in forming normaliser
 end;

 switch lower(OPT.alg.type)          % Calculate parameter update direction
  case 'lms',        direction = OPT.alg.mu*phi;
  case {'rls','ct'}, direction = Pphi/(OPT.alg.lambda+denom);   
  case 'kf',         direction = Pphi/(OPT.alg.R+denom);       
 end;
 
 th = th + direction*e;              % Update parameter estimate
 
 switch lower(OPT.alg.type)          % Update parameter covariance matrix
  case {'rls','ct'}, P = (P - Pphi*(Pphi')/(OPT.alg.lambda+denom))/OPT.alg.lambda;
  case {'kf'},       P = P - P*phi*(Pphi')/(OPT.alg.R+denom) + OPT.alg.Q;     
 end;
 if strcmp(lower(OPT.alg.type),'ct') P = (P/trace(P))*OPT.alg.P; end;

 % Save the time history of the parameter evolution.
 G.th_hist(t,:) = th'; 

 % Save the prediction error history 
 G.pe(t) = e;
end;

G.th = th;  % Output the final estimate.
G.P = P;    % Final Covariance Matrix    










