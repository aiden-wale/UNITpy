% This function tries to compute the maximum likelihood estimate of a
% state-space model from frequency domain measurements using the
% Expectation-Maximisation approach. This function is not meant to be
% called directly.  Rather, it is intended to be called by using the
% est.m algorithm, with OPT.alg='em' specified.
%
%     Usage is:  G = fem(Z,M,OPT);
%
%     where
%
%       Z:        Frequency-Response-Function (FRF) data. Z.y(i,j,k) holds
%                 the i'th output, j'th input FRF data for the k'th
%                 frequency point. The frequency points are stored in Z.w.
% M.ss.A,B,C,D:   Initial state space model structure guess.  These *cannot*
%                 be specified simply as integer orders, they must be an
%                 actual full state space system parameterisation.
%   M.ss.Q,R:     Covariance matrices for states and outputs, respectively.
%       M.op:     set to 'q' for shift and 's' for Laplace.  Default = 'q'.
%       M.T:      sampling period (ignored for q operator case).  Default=1
%       M.w:      vector of frequencies at which to calculate frequency
%                 response of estimated model.  Specify in real frequency,
%                 not normalised.  Default is 3 decades up to folding freq.
%      OPT:       Data structure which defines options for the estimation
%                 algorithm as follows:
%      OPT.miter: Maximum number of iterations in search for minimum.
%      OPT.dsp:   Control of output to screen 0=>quiet,1=>verbose.
%  OPT.stoptol:   Stopping criteria tolerance (1e-3 by default).
%      G:         Data structure which specifies the estimated model as
%                 follows:
%  G.ss.A,B,C,D:  [A,B,C,D] matrices defining estimated state space model.
%      G.G:       Matrix of frequency responses.  If the system has multiple
%                 inputs and multpile outputs, then this matrix is 3
%                 dimensional, with one `page' per output, and the i'th
%                 column of the j'th page (ie G.G(i,j,:)) being the
%                 frequency response from input i to ouput j.
%
%
%    Written by Brett Ninness, School of Elec. Eng and Comp. Sci.
%               Adrian Wills   University of Newcastle
%                              Australia.

% Copyright (C) Brett Ninness.

function g = fem(z,mm,opt);

%  Flag to we can make this same function work in both MATLAB or OCTAVE 
oct = exist('OCTAVE_VERSION');

% Set some options if they have not been specified
if ~isfield(opt,'miter')   opt.miter   = 300;        end; % default number of total iterations
if ~isfield(opt,'dsp')     opt.dsp     = 0;          end; % display iteration info (dsp = 1 to display)
if ~isfield(opt,'delta')   opt.delta   = 1e3;        end; % default trust region radius
if ~isfield(opt,'bismax')  opt.bismax  = 50;         end; % maximum number of trust region radius bisections
if ~isfield(opt,'stoptol') opt.stoptol = 1e-3;       end; % Newton decrement stop tolerance
if ~isfield(opt,'optit')   opt.optit   = 100;        end; % default number of optimisation iterations
if ~isfield(opt,'smeth')   opt.smeth   = 'bilin';    end; % default method for dealing with continuous models (can also be direct = use j*w)
if ~isfield(opt,'dspbode') opt.dspbode = 0;          end; % display bode plot at each iteration
if ~isfield(opt,'use_old') opt.use_old = 1;          end; % Old or new version of M step

% Detect if gui is running
gui = 0; guih = [];
if isfield(opt,'gui'),
 if ~isempty(opt.gui)
  gui  = 1;         %GUI is running
  guih = opt.gui;   %GUI handle
 end
end

% Extract data
z = startZ(z);
y = z.y;
w = z.w(:)';
N = max(size(y));

% Perform the pre-warp if continuous
if strcmp(mm.op,'s'),
 switch opt.smeth
  case 'bilin',
   T      = 2*pi/max(w);
   Tq     = 1;
   w      = 2*atan(w*T/2);
   mm.op  = 'q';
   opsave = 's';
   mm.T   = 1;
   ejw = exp(j*Tq*w); ejw = ejw(:).';

  otherwise % use a direct approach
   opsave='s';
   ejw = j*w; ejw = ejw(:).';
 end
else
 Tq = mm.T;
 opsave = 'q';
 ejw = exp(j*Tq*w); ejw = ejw(:).';
end

% Set initial values
del1 = opt.delta;
del2 = opt.delta;
a = mm.ss.A;
b = mm.ss.B;
c = mm.ss.C;
d = mm.ss.D;
n = size(a,1);
m = size(b,2);
p = size(c,1);
nnm=n*(n+m); nn=n*n;
if ~isfield(mm.ss,'Q'),
 Q = 10*eye(n);
else
 Q = mm.ss.Q;
end
if ~isfield(mm.ss,'R'),
 R = eye(p);
else
 R = mm.ss.R;
end

% Make input equal to the identity for each w(k) if not supplied by user
if ~isfield(z,'u'),
 z.u = zeros(m,m,N);
 for k=1:N, z.u(:,:,k) = eye(m); end
elseif isempty(z.u),
 z.u = zeros(m,m,N);
 for k=1:N, z.u(:,:,k) = eye(m); end
end
u = z.u;

% Modify system matrices for prewarp
if strcmp(opsave,'s') && strcmp(opt.smeth,'bilin'),
 a=pinv(eye(n)-(T/2)*a)*((T/2)*a+eye(n));
 b=(sqrt(T)/2)*(a*b+b);
 c=(sqrt(T)/2)*(c*a+c);
 d=d+c*inv(eye(n)+a)*b;
 udisp('Using Bilinear transform to handle continuous domain data.',gui,guih)
end

% Display iteration info if asked to
if opt.dsp,
 udisp('------------------------------------------------------------------------------------------',gui,guih);
 udisp(sprintf('%8s%20s%20s%10s%10s%20s','Iter #','Log-Likelihood','Pred. Error','Op. Its','No. Bis.','Newton-Decrement'),gui,guih);
 udisp('------------------------------------------------------------------------------------------',gui,guih);
end

if opt.dspbode,
 sysG = frd(z.y,w);
 if mm.op=='q'
  bode(sysG,ss(a,b,c,d,mm.T),w); drawnow;
 elseif mm.op=='s'
  bode(sysG,ss(a,b,c,d),w); drawnow;
 end
end

% start main EM loop
final_print = 0; %=1 prints final iteration information
numit       = 0; %number of iterations taken in minloop
kk          = 0; %number of bisections return from minloop
gnorm       = 0; %gradient norm returned from minimisation step
count       = 0; %counter for main loop
while count<opt.miter,
 count=count+1; %update main loop counter

 %keyboard
 
 % Call E-STEP
 %try
  [Phi,Psi,Sigma,LL(count),MSE(count)] = fem_estep(N,n,m,p,y,u,ejw,a,b,c,d,Q,R,mm.op,w);
 %catch
%   if opt.dsp, udisp('Could not complete E-Step',gui,guih); end
%   final_print=1;
%   break;
%  end

 % Check to see if we have gone backwards and quit if necessary
 if count>1 & LL(count)<LL(count-1),
  % Recall the last system before the decrease
  a=as; b=bs; c=ccs; d=ds; Q=qs; R=rs;
  if opt.dsp, udisp('Decrease in Log-Likelihood',gui,guih); end
  final_print=1;
  break;
 end

 % Check to see if we can stop because MSE has gone below user defined value
 if isfield(opt,'sv'),
  if MSE(count)<opt.sv,
   if opt.dsp, udisp('MSE less than OPT.sv',gui,guih); end
   final_print=1;
   break;
  end
 end

 % Save the system matrices in case we decrease the likelihood
 as=a; bs=b; ccs=c; ds=d; qs=Q; rs=R;
    
 % Reset delta's for each pass through minimisation loop
 del1=opt.delta; del2=opt.delta;

 % Call M-STEP
 %try
  [a,b,c,d,Q,R,del1,del2,kk,numit,gnorm] = fem_mstep(N,n,m,p,a,b,c,d,ejw,Phi,Psi,Sigma,opt,del1,del2,count);
%  catch
%   if opt.dsp, udisp('Could not complete M-Step',gui,guih); end
%   final_print=1;
%   break;
%  end

 % Provide iteration update if required
 if opt.dsp,
  udisp(sprintf('%8i%20.5e%20.5e%10i%10i%20.5e',count,LL(count),MSE(count),numit,kk,gnorm),gui,guih)
 end

 % Display bode plot if requested
 if opt.dspbode,
  if mm.op=='q'
   bode(sysG,ss(a,b,c,d,mm.T),w); drawnow;
  elseif mm.op=='s'
   bode(sysG,ss(a,b,c,d),w); drawnow;
  end
 end
 
end

% Provide iteration update if required
if opt.dsp & final_print,
 udisp(sprintf('%8i%20.5e%20.5e%10i%10i%20.5e',count,LL(count),MSE(count),numit,kk,gnorm),gui,guih)
end

if opt.dsp,
 udisp('------------------------------------------------------------------------------------------',gui,guih);
end

% Now put new things into g
if strcmp(opsave,'s') & strcmp(opt.smeth,'bilin'),
 aa   = inv(eye(size(a))+a);
 a    = (2/T) * aa * (a-eye(size(a)));
 d    = d - c*aa*b;
 b    = (2/sqrt(T))*aa*b;
 c    = (2/sqrt(T))*c*aa;
 mm.T = 0;
end
mm.op   = opsave;
g       = mm;
g.ss.A  = a;
g.ss.B  = b;
g.ss.C  = c;
g.ss.D  = d;
g.ss.Q  = Q;
g.ss.R  = R;
g.ss.K  = [];
g.ss.F  = [];
g.ss.G  = [];
g.LL    = LL;
g.MSE   = MSE;
try,
 if isfield(mm,'T'), g.T = mm.T; elseif strcmp(mm.op,'q'), g.T = 1; else g.T=0; end
 if isfield(mm,'w'), g.w = mm.w; else g.w = z.w; end
 g.op    = mm.op;
 g.delay = zeros(m,1);
 g.type  = 'ss';
catch,

end
%g = addtf(g);

g = sstotf(g);
g.var = 1;

% Add legend for prospective plotting
g.disp.legend=['Estimated ',g.type,' model'];
g.alg='EM'; % Record that block solution was used


% Finally make sure we give a Matlab system form
if ~oct
 if strcmp(mm.op,'s'),
  g.ss.sys=ss(g.ss.A,g.ss.B,g.ss.C,g.ss.D);
 else,
  g.ss.sys=ss(g.ss.A,g.ss.B,g.ss.C,g.ss.D,mm.T);
 end
end;

return;
%--------------------------------------------------------------------------
% END OF MAIN ROUTINE
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% START OF E-STEP ROUTINE
%--------------------------------------------------------------------------

function [Phi,Psi,Sigma,LL,MSE] = fem_estep(N,n,m,p,y,u,ejw,a,b,c,d,Q,R,op,w)
if strcmp(op,'q'), discrete=1;   else discrete=0;   end
[ve,de]=eig(a); de=diag(de); vei=ve\eye(size(ve));
[PP,Pssum,Psejw,Psej2,MSE,LL]=fsmooth(ve,de,vei,b,c,d,Q,R,w,ejw,y,u,discrete);
MSE=MSE/N; PP=PP*PP';
Phi=PP(1:n+p,1:n+p); Phi(1:n,1:n)=Phi(1:n,1:n)+Psej2;
Psi=PP(1:n+p,n+p+1:end); Psi(1:n,1:n)=Psi(1:n,1:n)+Psejw;
Sigma=PP(n+p+1:end,n+p+1:end); Sigma(1:n,1:n)=Sigma(1:n,1:n)+Pssum;
%--------------------------------------------------------------------------
% END OF E_STEP
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% START OF M-STEP ROUTINE
%--------------------------------------------------------------------------

function [a,b,c,d,Q,R,del1,del2,kk,numit,gnorm] = fem_mstep(N,n,m,p,a,b,c,d,ejw,Phi,Psi,Sigma,opt,del1,del2,itn);

if itn == 0,
 cdi = real(Psi(n+1:end,:))/real(Sigma);
 abi = real(Psi(1:n,:))/real(Sigma);
else
 cdi = [c d];
 abi = [a b];
end

% Maximise Q function with respect to A,B,C,D system matrices
if opt.use_old,
 [cd,del2,kk,numit,gnorm] = argmin('costcd',N,n,m,p,cdi,ejw,Phi(n+1:end,n+1:end),Psi(n+1:end,:),Sigma,opt,del2);
 [ab,del1,kk,numit,gnorm] = argmin('costab',N,n,m,p,abi,ejw,Phi(1:n,1:n),Psi(1:n,:),Sigma,opt,del1);
else
 [cd,numit,gnorm] = mincd(N,n,m,p,ejw,cdi,Phi(n+1:end,n+1:end),Psi(n+1:end,:),Sigma);
 [ab,numit,gnorm] = minab(N,n,m,p,ejw,abi,Phi(1:n,1:n),Psi(1:n,:),Sigma);
 kk   = 1;
 del1 = 1;
 del2 = 1;
end

% If no iteration, then make sure values are sensible
if opt.optit==0, kk=0; gnorm=0; end

% Combine the A B C and D matrices into G = Gamma matrix
G=[ab;cd];

% Form the covariance matrix Pi
CP = rchol([Phi -Psi;-Psi' Sigma]);
PP = CP*[eye(n+p);G'];
Pi = (PP'*PP)/N;

% Extract new system matrices and covariance matrices
a = G(1:n,1:n);
b = G(1:n,n+1:end);
c = G(n+1:end,1:n);
d = G(n+1:end,n+1:end);
Q = Pi(1:n,1:n);
R = Pi(n+1:end,n+1:end);

%----------------------------------------------------------------
% END OF M-STEP
%----------------------------------------------------------------

function [G,i,gnorm] = minab(N,n,m,p,ejw,G,Lambda,Omega,Pi);

[CC,pC]=chol([Lambda -Omega;-Omega' Pi]);
whilecnt = 0; alpC = 100*eps;
while pC~=0 & whilecnt<100,
    whilecnt = whilecnt + 1;
    [CC,pC]=chol([Lambda -Omega;-Omega' Pi] + alpC*eye(size(Lambda,1)+size(Pi,1)));
    alpC = 2*alpC;
end

for i=1:10,
    cc=CC*[eye(size(Lambda));G']; 
    R=(cc'*cc)/N; 
    Ri=pinv(R); 
    Ri=(Ri+Ri')/2; 
    gnorm = norm(vec(real(Ri*(G*Pi-Omega))));
    if gnorm > 1e-4, 
        HHH=(kron(real(Pi),real(Ri))-kron(imag(Pi)',imag(Ri)));
        gg=vec(real(Ri)*real(Omega)-imag(Ri)*imag(Omega)); 
        G=reshape(HHH\gg,n,n+m);
    else
        break;
    end
end

function [G,i,gnorm] = mincd(N,n,m,p,ejw,G,Phi,Psi,Sigma);

CC=chol([Phi -Psi;-Psi' Sigma] + 100*eps*eye(size(Phi,1)+size(Sigma,1)));

for i=1:10,
    cc=CC*[eye(size(Phi));G']; 
    R=(cc'*cc)/N; 
    Ri=pinv(R); 
    Ri=(Ri+Ri')/2; 
    gnorm = norm(vec(real(Ri*(G*Sigma-Psi))));
    if gnorm > 1e-4, 
        HHH=(kron(real(Sigma),real(Ri))-kron(imag(Sigma)',imag(Ri))); 
        gg=vec(real(Ri)*real(Psi)-imag(Ri)*imag(Psi)); 
        G=reshape(HHH\gg,p,n+m);
    else
        break;
    end
end


%----------------------------------------------------------------
% START of cost ab routine
%----------------------------------------------------------------
function [c,g,H1,H2] = costab(N,n,m,p,ejw,G,Lambda,Omega,Pi);

%Depending on whether or not we need derivatives, then compute eigenvectors
if nargout>1,
 [ve,de]=eig(G(:,1:n));
 de=diag(de);
else
 [de]=eig(G(:,1:n));
end

%Compute cost
c=cmpcost(ejw,de);

MM = [Lambda -Omega;-Omega' Pi];
MM = (MM+MM')/2;

[CC,pC]=chol(MM);
whilecnt = 0; alpC = 100*eps;
while pC~=0 & whilecnt<100,
    whilecnt = whilecnt + 1;
    [CC,pC]=chol(MM + alpC*eye(size(Lambda,1)+size(Pi,1)));
    alpC = 2*alpC;
end
%keyboard

CC=CC*[eye(size(Lambda));G'];
Q=(CC'*CC)/N;
c=real(c+N*log(det(Q)));

%Compute gradient
if nargout>1,
 nnm     = n*(n+m); 
 nn      = n*n;

 %Compute gradient vector
 vei     = ve\speye(size(ve)); 
 v1      = kronaw(ve,vei.');
 [J2,H2] = cmpgh(ejw,de); 
 J2      = v1*J2;
 g1      = 2*real(vec(pinv(Q)*[G*Pi-Omega]));
 g2      = 2*real(J2);
 g       = g1;
 g(1:nn) = g(1:nn)+2*real(J2);
end

if nargout>2,
 %Hack
 Qi=Q\speye(n);
 Qi=(Qi+Qi')/2;
 H1=2*real(kronaw(Pi.',Qi));

 %Compute Hessian matrix
 v2=kronaw(vei,ve.');
 idx=vec(reshape([1:nnm],n+m,n)');
 t1=kronaw(eye(n),[G*Pi-Omega]);
 t2=kronaw([G*Pi.'-conj(Omega)],eye(n));
 J=(t1(:,idx)+t2)/N;
 QikQi=kronaw(Qi,Qi.');

 scal=H2(:); scal=scal(:,ones(1,nn));
 HH=2*real(v1*(scal.*v2));
 idx=vec(reshape([1:nn],n,n)');
 HH=HH(:,idx);

 H2=-N*real(J.'*QikQi*conj(J));
 H2(1:nn,1:nn)=H2(1:nn,1:nn)+HH;
end
%----------------------------------------------------------------
% END of cost ab routine
%----------------------------------------------------------------

%----------------------------------------------------------------
% START of cost cd routine
%----------------------------------------------------------------
function [c,g,H1,H2] = costcd(N,n,m,p,ejw,G,Phi,Psi,Sigma);

%Compute cost
c=N*real(log(det(Phi-G*Psi' - Psi*G' + G*Sigma*G')));

if nargout>1,
 CC=rchol([Phi -Psi;-Psi' Sigma] + 100*eps*eye(size(Phi,1)+size(Sigma,1)));
 CC=CC*[eye(size(Phi));G'];
 R=(CC'*CC)/N;
 Ri=R\speye(p);
 Ri=(Ri+Ri')/2;
 gsp=G*Sigma-Psi;
 g=2*vec(real(Ri*gsp));
end

if nargout>2,
 pnm=p*(n+m); idx=vec(reshape([1:pnm],n+m,p)');
 t1=kronaw(eye(p),gsp);
 t2=kronaw(conj(gsp),eye(p));
 J=(t1(:,idx)+t2)/N; RikRi=kronaw(Ri,Ri.');
 H1=2*real(kronaw(Sigma.',Ri));
 H2=-N*real(J.'*RikRi*conj(J));
end


%----------------------------------------------------------------
% START of argmin routine
%----------------------------------------------------------------
function [G,del,kk,numit,gnorm] = argmin(costfn,N,n,m,p,G,ejw,Phi,Psi,Sigma,opt,del);

%Get size of parameter matrix
[nr,nc] = size(G);

%Stopping tolerance on Newton decrement
stoptol = opt.stoptol;

%Loop through Newton search opt.optit number of times
Ith=eye(length(vec(G))); H=Ith;
for numit=1:opt.optit,

 %Form gradient and Hessian
 [cost_old,g,H1,H2] = feval(costfn,N,n,m,p,ejw,G,Phi,Psi,Sigma);
 H=H1+H2;
 H=real(H+H')/2;


 %Use trust region method.
 S=eig(H); if ~isreal(S), udisp('Hessian has complex eigenvalues',gui,guih); end
 lam_min=min(0,min(real(S)));

 if lam_min < 0, lam=-2*lam_min; else, lam_min=0; lam=0; end
 eta=0.01;

 for kk=1:100,
  %Get search direction
  [cH,pp]=chol(H+lam*eye(size(H)));
  while pp~=0,
   lam=2*(lam+1e-3);
   [cH,pp]=chol(H+(lam+eps)*eye(size(H)));
  end
  dp=-(cH\(cH'\g));

  %Check that search direction satisfies ||dp|| <= del,
  %if not then find Lagrange multiplier for constrained
  %problem.
  i=0;
  if norm(dp)>del,
   hitdel=1;
   for i=1:100,
    q=cH'\dp;
    dlam=(norm(dp)/norm(q))^2*((norm(dp)-del)/del);
    lam=max(-lam_min,lam+dlam);
    [cH,pp]=chol(H+(lam+eps)*eye(size(H)));
    while pp~=0,
     lam=2*(lam+1e-1);
     [cH,pp]=chol(H+(lam+eps)*eye(size(H)));
    end
    dp=-cH\(cH'\g);
    if abs(norm(dp)-del) < 1e-3, break; end
   end
  else
   hitdel=0;
  end


  %Check stopping criteria
  gnorm=-g'*dp/length(g);
  if kk>5 && gnorm < stoptol, break; end

  %Get new cost
  Gp=reshape(dp,nr,nc); Gn=G+Gp;
  cost_new = feval(costfn,N,n,m,p,ejw,Gn,Phi,Psi,Sigma);

  %Compute actual/predicted cost reduction
  rho = (cost_old - cost_new)/(eps-g'*dp-0.5*(dp'*(H*dp + lam*dp)));

  %Adaptively change trust region according to local
  %performance.
  if rho < 0.25,
   del=0.25*norm(dp);
  elseif rho > 0.75 & hitdel,
   del=min(2*del,1e16);
  end

  %If we reduced the cost sufficiently then stop.
  if rho > eta, G=Gn; break; end
 end

 %Break if we have a sufficiently small Newton decrement
 %gnorm=abs(dp'*g)/length(g);
 if gnorm < stoptol, break; end
 if norm(dp)/norm(G(:)) < 1000*eps, break; end
end
%----------------------------------------------------------------
% END of argmin routine
%----------------------------------------------------------------
