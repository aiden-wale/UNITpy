%  Function used as part of DDLC method for estimation of state-space
%  models. Function converts from a minimally parametrised vector beta to
%  a fully parametrised vector theta via the equation
% 
%    theta = theta_0 + Qp*beta. 
%
%  After theta has been formed as above, a matrix Q is constructed whose
%  column range space specifies directions (in theta space) that correspond
%  (locally) to similarity transformations of the state-space system. The
%  column null-space matrix of Q is calculated and stored in Qp, hence the
%  columns of Qp describe directions which are orthogonal to directions of
%  similarity transformation.
%
%  The state-space system is assumed to be of the (possibly bilinear) form
%
%  x(k+1)  =  A*x(k) + F*(u(k) kron x(k)) + B*u(k)
%
%    y(k)  =  C*x(k) + G*(u(k) kron x(k)) + D*u(k)
%
%  where kron denotes the Kronecker-Tensor product. With this in mind,
%  the column vector theta is parametrised according to
%
%           +        +
%           | vec(A) |
%           | vec(B) |
%           | vec(C) |
%  theta  = | vec(D) |
%           | vec(K) |
%           | vec(F) |
%           | vec(G) |
%           +        +
%
%  Usage is: 
%
%  [G,beta_new] = ddlc(beta,M);
%
%  Where
%
%  M        = Model structure definition in MATLAB structure.
%  M.ss.A   = System A matrix (same for B,C,D,K,F,G). Note that D,K,F
%             and G may be empty if desired.
%  M.th0    = Initial theta vector (theta_0) as described above.
%  M.Qp     = Matrix as described above (dimension is n_theta by n_beta).
%  beta     = Minimal parametrisation (same number of elements as the
%             column dimension of M.Qp).
%
%  G        = Copy of M but with new calculated G.Q and G.Qp corresponding
%             to new theta. G.ss.(A,B,C,D,K,F,G) are also update according to
%             theta. Note that G.ss.F and G.ss.G may be set to empty matrices if
%             a purely linear system is of interest;
%  G.th0    = theta as calculated above.
% G.Q,G.Qp  = Matrices with (respectively) columns spanning the tangent space
%             of systems at the point G.th0, and the orthogonal complement
%             of the tangent space of systems at G.th0;
%  beta_new = a column vector of zeros so as to ensure that
%             theta = G.th0 + G.Qp*beta_new. 
%
%
% written by  Brett Ninness, School of EE & CS
%             Adrian Wills   University of Newcastle
%                            Australia.

% Copyright (C) Brett Ninness

function [m]=ddlc(th,M);

%Check for bilinear parts and set to empty if missing
if ~isfield(M.ss,'F'),
 M.ss.F=[];
end
if ~isfield(M.ss,'G'),
 M.ss.G=[];
end
M = theta2m(th,M,1);

% Return new system in both full ss and m.th0 + m.Qp*theta (zeta = 0) form
m=M;

% Work out dimensions of state, input, output, etc.
nx = size(M.ss.A,1); 
nu = size(M.ss.B,2); 
ny = size(M.ss.C,1); 
In = eye(nx);
if ~isempty(M.ss.D) 
 zz=zeros(nu*ny,nx^2); 
else
 zz=[]; 
end

% Columns of Q span tangent space (in theta domain) of equivalent systems
Q=[kron(M.ss.A',In)-kron(In,M.ss.A);kron(M.ss.B',In);-kron(In,M.ss.C);zz;kron(M.ss.K',In)];
if ~isempty(M.ss.F),for i=1:nu,Q=[Q;kron(M.ss.F(:,(i-1)*nx+1:i*nx)',In)-kron(In,M.ss.F(:,(i-1)*nx+1:i*nx))];end;end
if ~isempty(M.ss.G),for i=1:nu,Q=[Q;-kron(In,M.ss.G(:,(i-1)*nx+1:i*nx))];end;end
m.Q = Q;

% Find basis for space orthogonal to this
[P,R]=qr(Q); m.Qp=P(:,size(Q,2)+1:end);
