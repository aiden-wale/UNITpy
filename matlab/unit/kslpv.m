%
% KS. This function is the square-root implementation of the Kalman
% Smoother as per
%
% [1] Gibson and Ninness, "RobustMaximum-Likelihood Estimation of Multivariable
%     Dynamic Systems", Automatica, , 41(10):1667?1682, October 2005.
%
% *NOTE*: Covariance information is returned in SQUAREROOT form
%
% Model is assumed to be of the form,
%
% x(t+1) = A(t)x(t) + B(t)u(t) + w(t),   [w(t)]    (    [ Q(t)    S(t) ] )
%                                        [    ] ~ N( 0, [              ] )
%   y(t) = C(t)x(t) + D(t)u(t) + v(t),   [v(t)]    (    [ S^T(t)  R(t) ] )
%
% A call to this function should look like
%
%       G = ks(Z,M,OPT);
%
% Where
%
%            Z:  Input-Output data in one of two forms.  The standard form
%                is for it to be a record with elements Z.y and Z.u, each
%                of which are matrices with number of rows equal to the
%                number of data samples, and number of columns equal (respectively)
%                to the number of outputs and the number of inputs.  On
%                the other hand, Z can be a matrix of the form Z = [y,u]
%                where it is assumed that y is a column vector of output
%                measurements and u is a matrix whose columns are the
%                input measurements; in this latter MISO models are
%                being considered.
%
%             M: Data structure which defines the above model:
%  M.ss.A,B,C,D: Possibly time-varying system matrices; for each matrix, it is
%                assumed that the time index is the third dimension, e.g.
%                A(:,:,t) is the state transition matrix at time t. If the
%                third dimension is equal to one, then time-invariant
%                matrices are assumed for that case.
%
%    M.ss.Q,S,R: Possibly time-varying noise covariance matrices; for each
%                matrix, it is assumed that the time index is the third
%                dimension. If the third dimension is equal to one, then
%                time-invariant matrices are assumed for that case.
%
%    M.ss.X1,P1: Initial state mean (X1) and its covariance matrix (P1),
%                respectively.
%
%             G: returned structure with the following fields
%          G.xp: predicted states, i.e. E[x(t) | y_1,..,y_{t-1}]
%          G.xf: filtered states, i.e.  E[x(t) | y_1,..,y_{t}]
%          G.xs: smoothed states, i.e.  E[x(t) | y_1,..,y_{N}]
%
%          G.Pp: *SQUAREROOT* of predicted state covariance matrix, i.e.
%                  Pp(:,:,t) = E{x(t)*x(t)' |  y_1,..,y_{t-1}}
%          G.Pf: *SQUAREROOT* of filtered state covariance matrix
%                  Pf(:,:,t) = E{x(t)*x(t)' |  y_1,..,y_{t}}
%          G.Ps: *SQUAREROOT* of smoothed state covariance matrix, i.e.
%                  Ps(:,:,t) = E{x(t)*x(t)' |  y_1,..,y_{N}}
%          G.Ms: cross covariance between x(t+1) and x(t), i.e.
%                  Ms(:,:,t) = E{x(t+1)*x(t)' | y_1,..,y_{N}}
%          
%          G.yp: predicted output estimate
%          G.yf: filtered output estimate
%          G.ys: smoothed output estimate
%          G.pe: prediction error
%          G.fe: filter error
%          G.se: smoother error
%
%          G.LL: negative log-likelihood
%
%   written by Adrian Wills:  School of EE & CS
%              Brett Ninness  University of Newcastle
%        	   Stuart Gibson  Australia.
%               Last Revised  20/03/2009.
%

% Copyright (C) 1999-2009 Brett Ninness.

function G = ks(Z,M,OPT)


%Now get default options
if nargin<3, OPT = []; end
OPT = startOPT(OPT);

%Extract data and sizes
[y,u,p,m,N] = Z2data(Z);
y = y';

%If we have no measurements then we cant filter or smooth
if p<1,
	error('According to the data, there are no outputs. Nothing to do.')
end

% Figure out what parts of model are specified and set the rest to defaults
if ~exist('M')
	error('Need to specify initial model structure M!');
elseif isfield(M,'ss')
	if ~isfield(M.ss,'A'),   error('Need to specify M.ss.A!');              end
	%Get state dimension
	n = size(M.ss.A,1);
	if ~isfield(M.ss,'C'),   error('Need to specify M.ss.C!');              end
	if ~isfield(M.ss,'B'),   M.ss.B  = zeros(n,m);                          end
	if ~isfield(M.ss,'D'),   M.ss.D  = zeros(p,m);                          end
	if ~isfield(M.ss,'R'),   M.ss.R  = eye(p,p);                            end
	if ~isfield(M.ss,'Q'),   M.ss.Q  = 0.001*eye(n,n);                      end
	if ~isfield(M.ss,'S'),   M.ss.S  = zeros(n,p);                          end
	if ~isfield(M.ss,'X1'),  M.ss.X1 = zeros(n,1);                          end
	if ~isfield(M.ss,'P1'),  M.ss.P1 = eye(n,n);                            end
	if ~isfield(M,'op'),     M.op    = 'q';                                 end
	if ~isfield(M,'T'),      M.T     = 1;                                   end
	if ~isfield(M,'delay'),  M.delay = zeros(m,1);                          end
	if ~isfield(M,'type'),   M.type  = 'ss';                                end
else
	error('Need to specify model in M.ss fields');
end

%Check for input or not and adjust everything else to suit
if m>0,
	% Include delays specified in model structure on inputs
	for r=1:m
		u(:,r) = [zeros(M.delay(r),1);u(1:N-M.delay(r),r)];
	end
	u = u';
else
	u = zeros(0,N);
end

% Run some checks for bilinear systems
if m<1 && ~isempty(strfind(M.type,'bilin')), % If there is no input, then the system is not bilinear between input and state, which is all we support
	M.type = 'ss'; 
end

%Make sure we set M.type to bilinear if it is bilin
if strfind(M.type,'bilin'),
	M.type = 'bilinear';
end

% Set OPT.smoothing so that smoothing is done
OPT.smoothing = 1;

% Now call the Kalman Filter/Smoother routine
Z.y = y;
Z.u = u;
G   = rksqrtv_lpv(Z,M,OPT);

% For the purposes of backward compatibility
G.ss.X  = G.ss.xs;
G.ss.Xt = G.ss.xs;
G.ss.Xs = G.ss.xs;
G.ss.Xf = G.ss.xf;
G.ss.Xp = G.ss.xp;
