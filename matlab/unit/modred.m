% MODRED - Function to find a reduced order model for a LTI system
% in state space form by using Kung's method.
%
% Usage is:
%
% [a,b,c,d]=modred(A,B,C,D,ord)
%
% Where
%
% A,B,C,D = System matrices defining LTI state space system
%     ord = Order required for reduced order model.  If not 
%           specified then will be chosen automatically via 
%           thresholding on singular values of impulse 
%           response Hankel matrix
%
% a,b,c,d = System matrices defining reduced order LTI state
%           space system model.
%
% written by Brett Ninness, School of EE & CS
%            Adrian Wills   University of Newcastle
%      		                Australia.

% Copyright (C) Brett Ninness

function [a,b,c,d]=modred(a,b,c,d,r);

%Extract sizes of things
m=size(b,2); p=size(c,1); n=size(a,1);

if nargin<5, r=0; end

if r>n, error('New order is greater than old one.'); end

%Return early if state dimension==1
if n<2, return; end

%Check to see if we have an input or not
if m==0 | b==0, input=0; else, input=1; end

%Build (ex)tended observability and controlability matrices
O=obsv(a,c); if input C=obsv(a',b')'; else, C=1; end;

%Use Kung's method to extract system matrices via balanced realization
[U,S,V]=svd(O*C); s=diag(S);

%Try and determine rank from singular values
if r==0, r=sum(s>s(1)*100*eps); end

%Make new Observability and Controlability matrices
O=U(:,1:r)*diag(sqrt(s(1:r))); 
if input, C=diag(sqrt(s(1:r)))*V(:,1:r)'; end

%Extract new system matrices (d stays the same)
a=O(1:end-p,:)\O(p+1:end,:); 
if input, b=C(:,1:m); end
c=O(1:p,:);