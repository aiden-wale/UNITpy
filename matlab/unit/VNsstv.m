function [cost,pe,grad,phi] = VNsstv(Z,theta,OPT,M,div)

y  = Z.y;           % Outputs
u  = Z.u;           % Inputs
t  = Z.t;           % Time stamps
d  = Z.d;           % Integration time(s)
nt = length(theta); % Number of parameters

% Rotate data matrices into correct form
if size(y,1)>size(y,2), y=y'; end
if size(u,1)>size(u,2), u=u'; end
[p,N] = size(y); 
N     = N-1;
m     = size(u,1);
if length(d)==1, d=d(ones(1,N)); end

% Set return variables to default values
cost = 0;
pe   = [];
grad = [];
phi  = [];

% Get state dimension
n    = size(M.ss.A,1);
deps = eps^(1/2);
x    = zeros(n,1);
P    = zeros(n);
pes  = zeros(Z.ny,Z.Ny);

% If we want derivatives and Hessian approx, then make some room
if div,
 phi  = zeros(N+N*p+N*p,nt);
 grad = zeros(nt,1);
 dP   = zeros(n,n,nt);
 dx   = zeros(n,nt);
end

% Finite difference gradient approximation method

findiff = 1;  %1 = forward difference,  2 = mid point

% Now run Kalman Filter
for k=1:N,
 [A,B,C,D,Q,S,R] = sample(M,theta,t(k+1)-t(k),d(k));
 pe              = y(:,k) - C*x - D*u(:,k);
 Lam             = C*P*C'+R;
 K               = (A*P*C'+S)/Lam;
 pes(:,k)        = pe;
 
 % Update cost
 detLam = det(Lam);
 if real(detLam)<=0, cost = inf; return; end
 Lpe    = Lam\pe;
 cost   = cost + pe'*Lpe + log(detLam);
 
 % Now compute derivatives if required
 if div,
  for i=1:nt,
   % Obtain derivatives of matrices via numerical differentiation
   dtheta    = zeros(nt,1);
   dtheta(i) = deps;
   [Ap,Bp,Cp,Dp,Qp,Sp,Rp] = sample(M,theta+dtheta,t(k+1)-t(k),d(k));
   if findiff==2
    [Am,Bm,Cm,Dm,Qm,Sm,Rm] = sample(M,theta-dtheta,t(k+1)-t(k),d(k));
   else
    Am=A; Bm=B; Cm=C; Dm=D; Qm=Q; Sm=S; Rm=R;
   end
   dA        = (Ap-Am)/(findiff*deps);
   dB        = (Bp-Bm)/(findiff*deps);
   dC        = (Cp-Cm)/(findiff*deps);
   dD        = (Dp-Dm)/(findiff*deps);
   dQ        = (Qp-Qm)/(findiff*deps);
   dS        = (Sp-Sm)/(findiff*deps);
   dR        = (Rp-Rm)/(findiff*deps);
   
   dL        = dC*P*C'+C*dP(:,:,i)*C' + C*P*dC' + dR;
   dK        = (dA*P*C' + A*dP(:,:,i)*C' + A*P*dC' + dS - K*dL)/Lam;
   dP(:,:,i) = dA*P*A' + A*dP(:,:,i)*A' + A*P*dA' + dQ - dK*Lam*K' - K*dL*K' - K*Lam*dK'; 
   
   de        = -dC*x - C*dx(:,i) - dD*u(:,k);
   dx(:,i)   =  dA*x + A*dx(:,i) + dB*u(:,k) + dK*pe + K*de;
   
   % Update gradient
   LdL       = Lam\dL;
   grad(i)   = grad(i) + trace(LdL) + 2*de'*Lpe - Lpe'*dL*Lpe;
   
   % Update phi
   phi(k,i)                         = trace(LdL);
   phi(N+(k-1)*p+1:N+k*p,i)         = 2*(chol(Lam)'\de);
   phi(N+N*p+(k-1)*p+1:N+N*p+k*p,i) = vec(LdL);
  end
 end
 
 % Update the state and covariance
 P = A*P*A' + Q - K*Lam*K';
 x = A*x + B*u(:,k) + K*pe;
end

% Normalise cost for data length
cost = cost/N;

% Do rotations for phi if divs asked for
if div,
 grad = grad/N;
 phi  = triu(qr(phi/sqrt(N),0)); 
 phi  = phi(1:nt,1:nt);
 pe   = phi'\grad(:);
else
 pe = pes';
end
