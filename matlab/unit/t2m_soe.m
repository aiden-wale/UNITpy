%  T2M_SOE: Function to allow for estimation of continuous time
%  transfer function structures by grey-box parametrization of
%  a state space structure - by moving to ss we know how to integrate
%  between possibly irregular spaced samples.
%
%  This function specifies the appropriate grey box mapping.
%
%  Usage is
%
%  G = t2m_soe(M,theta);
%
%  where
%
%  M     = Model structure specification.  It should have M.op='s' and 
%          M.type = 'oe' for this function to be getting called;
%
%  theta = Vector of parameters specifying an OE structure
%
%  G     = Structure with G.ss terms set to canonical state space model 
%          consistent with what is specified in theta.  Ie - OE model
%          in state space form.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness.

function M = t2m_soe(M,theta)

theta=theta(:); % Make sure theta is a column vector

% Extract current OE=num/den tf estimates from theta
num = theta(1:M.nB+1)'; den = [1,theta(M.nB+2:end)'];

% Express now as canonically parametrized ss model
[M.ss.A,M.ss.B,M.ss.C,M.ss.D] = tf2ss(num,den);

% Since it is an OE model, there is no state noise
nx = size(M.ss.A,1);
M.ss.Q  = zeros(nx,nx);
M.ss.K  = zeros(nx,1);
M.ss.S  = zeros(nx,M.ny);

% We are not going to estimate innovations variance
% in the gradient based search loop
M.ss.R  = zeros(M.ny,M.ny);

% We are not estimating initial state - but could very simply
% do so if we wanted.
M.ss.X1 = zeros(nx,1);