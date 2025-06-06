function G = kfonestep(A,B,C,D,Q,S,R,X1,P1,y,u)
%
% Square-root implementation of Kalman Filter/Smoother as per
%
% Gibson and Ninness, "RobustMaximum-Likelihood Estimation of Multivariable
% Dynamic Systems", Automatica, , 41(10):1667?1682, October 2005.
%
% Model is assumed to be of the form,
%
% x(t+1) = Ax(t) + Bu(t) + w(t),   [w(t)]    (    [ Q    S ] )
%                                  [    ] ~ N( 0, [        ] )
%   y(t) = Cx(t) + Du(t) + v(t),   [v(t)]    (    [ S^T  R ] )
%
% A call to this function should look like
%
%       G = kfonestep(A,B,C,D,Q,S,R,X1,P1,y,u);
%
% Where
%
%           y,u: are, respectively, column vectors of the output and input 
%                data for the current time t.
%
%       A,B,C,D: system matrices.
%
%         Q,S,R: covariance matrices.
%
%            X1: mean value of initial state
%            P1: *squareroot* of initial state covariance (upper triangular
%                squareroot, e.g. obtained via chol())
%
%             G: returned structure with the following fields
%
%          
%          G.xp: predicted states, i.e. E[x(t+1) | y(t)]
%          G.xf: filtered states, i.e.  E[x(t)   | y(t)]
%
%          G.Pp: *squareroot* of predicted state covariance matrix
%          G.Pf: *squareroot* of filtered state covariance matrix
%          
%          G.yp: predicted output estimate
%          G.yf: filtered output estimate
%          G.pe: prediction error
%          G.fe: filter error
%          G.K : Kalman Gain matrix
%          G.Ri: innovations covariance
%          G.LL: negative log-likelihood
%          G.PE: prediction error cost
%
%
%
%
%   written by Adrian Wills: School of EE & CS
%                            University of Newcastle
%        	                 Australia.
%               Last Revised 23/07/2009.
%

%Extract data and sizes
n = size(A,1);
p = size(y,1);
m = size(u,1);
N = 1;

%If we have no measurements then we can't filter
if p<1,
	error('According to the data, there are no outputs. Nothing to do.')
end

%Make some room for mean and covariance
G.xp = zeros(n,1);
G.xf = zeros(n,1);

G.Pp = zeros(n,n);  
G.Pf = zeros(n,n);
G.Ri = zeros(p,p);
G.K  = zeros(n,p);

G.yp = zeros(p,1);
G.yf = zeros(p,1);
G.pe = zeros(p,1);
G.fe = zeros(p,1);

%Make some room
G.LL  = 0;
G.PE  = 0;
R1    = zeros(n+p);
R2    = zeros(n+p);
R3    = zeros(2*n,n);
R4    = zeros(2*n,n);


% Robust computation of Q-S*inv(R)*S' that preserves symmetry and
% non-negativity
X = triu(rchol([R S'; S Q]));
R = X(1:p,1:p);
Q = X(p+1:p+n,p+1:p+n);
normS = norm(S(:));
if normS>0.0,
	SR1 = (X(1:p,p+1:p+n)')/(X(1:p,1:p)');
else
	SR1 = zeros(n,p);
end
if normS>0.0,
	A = A-SR1*C;
	B = B-SR1*D;
end
R1(1:p,1:p)         = X(1:p,1:p);
R1(p+1:end,1:p)     = P1*C';
R1(p+1:end,p+1:end) = P1;
R2                  = triu(qr(R1));
G.Ri                = R2(1:p,1:p);
G.K                 = R2(1:p,p+1:p+n)'/(R2(1:p,1:p)');
G.Pf                = R2(p+1:p+n,p+1:p+n);
G.yp                = C*X1 + D*u;
G.pe                = y - G.yp;
G.PE                = G.pe'*G.pe;
Riep                = rfbs(G.pe,R2(1:p,1:p)',1);
G.LL                = Riep(:)'*Riep(:) + 2*sum(log(abs(diag(R2(1:p,1:p)))));
G.xf                = X1 + R2(1:p,p+1:end)'*Riep;	
G.yf                = C*G.xf + D*u;
G.fe                = y - G.yf;
G.xp                = A*G.xf + B*u  + SR1*y;
R3(1:n,:)           = G.Pf*A';
R3(n+1:end,:)       = X(p+1:p+n,p+1:p+n);
R4                  = triu(qr(R3,0));
G.Pp                = R4(1:n,:);

return;

%--------------------------------------------------------------------------
%
%  AUXILIARY FUNCTIONS
%
%--------------------------------------------------------------------------
% Function to compute Cholesky factor robustly
function [A] = rchol(A)
A = triu(A); n = size(A,1); tol = n*eps;
if A(1,1) <= tol,
	A(1,1:n) = 0;
else
	A(1,1:n) = A(1,1:n)/sqrt(A(1,1));
end
for j=2:n,
	A(j,j:n) = A(j,j:n) - A(1:j-1,j)'*A(1:j-1,j:n);
	if A(j,j) <= tol,
		A(j,j:n) = 0;
	else
		A(j,j:n) = A(j,j:n)/sqrt(A(j,j));
	end
end


% Function that performs robust forward or backward substitution
function X = rfbs(B,A,uplow)
%uplow = 1 for lower triangular, 0 for upper
deps=100*eps; [n,m]=size(B); X=zeros(n,m);
if uplow
	for i=1:m,
		if abs(A(1,1))>deps,
			X(1,i)=B(1,i)/A(1,1);
		end
		for j=2:n,
			if abs(A(j,j))>deps,
				X(j,i)=(B(j,i)-A(j,1:j-1)*X(1:j-1,i))/A(j,j);
			end
		end
	end
else
	for i=1:m,
		if abs(A(n,n))>deps,
			X(n,i)=B(n,i)/A(n,n);
		end
		for j=n-1:-1:1,
			if abs(A(j,j))>deps,
				X(j,i)=(B(j,i)-A(j,j+1:n)*X(j+1:n,i))/A(j,j);
			end
		end
	end
end