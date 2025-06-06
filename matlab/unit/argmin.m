%  Function to find minimum of function that is defined by user supplied
%  m file.  Algorithm used is a Damped Newton direction line-search, possibly
%  switching to a Gradient line search if the Newton direction is unable to
%  decrease the cost on a particular iteration.
%
%  Usage is
%
%  [theta,cost_log,args] = argmin(Z,'FUN',theta0,OPT,ARGS)
%
%  where
%
%  Z           = data that cost defined by `FUN' depends upon.
%  FUN         = Handle for function that accepts Z and theta,OPT,ARGS
%                (in that order) and returns the cost, prediction error sequence,
%                gradient and Jacobian (in that order) at the point
%                theta.  If OPT.cmpgrad = 1, then gradient and Jacobian
%                will be computed using finite differences and only cost
%                and prediction error sequence are required.
%  theta0      = initial guess for minimising vector argument.
%  OPT.lmax    = maximum number of bisections in line search direction. Default = 52.
%  OPT.tol     = minimum norm on gradient vector before stopping search. Default = 1e-6.
%  OPT.miter   = maximum number of iterations in search for minimum.  Default = 100.
%  OPT.dsp     = 0 => quiet, 1 => verbose output.  Default = 0
%  OPT.mdec    = Minimum relative decrease of cost before search is
%                terminated.  Default = 1e-4.
%  OPT.subtol  = Relative sngular value tolerance, which determines how many singular
%                values (r) to include in the Hessian approximation according to the
%                rule r=max(i : s(i) > s(1)*OPT.subtol), where s(i) are singular
%                values in descending order (where s(1) is the largest). Default = 1e-6.
%  OPT.adapt   = If set to 1 (default) then the relative singular value tolerance
%                OPT.subtol is selected adaptively according to algorithm progress.
%                If set to 0 then the relative singular value tolerance OPT.subtol
%                is held at a constant value throughout the algorithm. Default = 1;
%  OPT.cmpgrd  = If set to 1 will numerically evaluate gradient and
%                Jacobian so that FUN does not have to supply them, 
%  OPT.ngt     = If set to 1 will numerically evaluate the gradient and
%                display it relative to the analytical gradient supplied
%                by FUN.  Handy for debugging.
%  OPT.saveit  = If set to 1 then the parameters values theta_k are saved
%                in ARGS.thetait as the k'th column upon return.
%  ARGS        = Structure which is passed to FUN and GRADFUN in case
%                auxiliary variables are required by them.
%
%  theta       = Terminal value of vector at end of minimisation search.
%                That is, theta is (hopefully) the minimising argument of
%                the specified cost function;
%  cost_log    = Vector recording history of how cost was gradually
%                decreased as iterative search proceeded;
%  args        = If optional ARGS are passed to cost function, which may
%                be modified by optional function fargs at end of each
%                line-search iteration, then their terminal value is
%                passed back out as args.
%
%   written by Brett Ninness  School of EE & CS
%              Adrian Wills   University of Newcastle
%        		              Australia.


% Copyright (C) Brett Ninness, Adrian Wills

function [theta,cost_log,ARGS] = argmin(Z,FUN,theta,OPT,ARGS,fargs)

% Figure out what has not been specified in OPT structure and set to
% defaults
OPT = startOPT(OPT);

%Detect if gui is running
gui = 0; guih = [];
if isfield(OPT,'gui'),
 if ~isempty(OPT.gui)
  gui  = 1;         %GUI is running
  guih = OPT.gui;   %GUI handle
 end
end

%make sure data is in correct format
Z=startZ(Z);

%Start the display of information if desired
if OPT.dsp, 
 udisp('Algorithm: Gradient Based Search',gui,guih);
 switch OPT.cost,
  case {'det','ml'}
   mincost = 'Maximum Likelihood';
  otherwise
   mincost = 'Mean Squared Error';
 end
 udisp(['Cost: ',mincost],gui,guih);
end
sglines='----------------------------------------------------------------------------';

% Do we want to save iterates or not?
if ~isfield(OPT,'saveit'),
 OPT.saveit = 0;
elseif OPT.saveit,
 theta_saved = [theta];
end

% determiine if we have a theta structure
if isfield(ARGS,'theta_struct'),
 theta_struct = ARGS.theta_struct;
else
 theta_struct = find(ones(size(theta)));
end

% Initialisation
count     = 1; 
cost_log  = []; 
svtol     = OPT.subtol; 
delta     = OPT.delta; 
delta_rgn = 0;
trust_1st = 1;
Ith       = eye(length(theta_struct)); 
H         = Ith; 

hess_supplied = 0;
%If we already have a good inv Hessian then use it
if ~isempty(ARGS),
 if isfield(ARGS,'invHess'),
  H = ARGS.invHess;
  hess_supplied = 1;
 end
end


%If the user wants things displayed then we should popup the stop
%estimation GUI so that they can stop the performance...
if ~gui & OPT.dsp,
 try
  stopEstGui=stopEstimation();
  stopEstGui.setVisible(true);
  stopEstGuiRunning = 1;
 catch
  stopEstGuiRunning = 0;
 end
else
 stopEstGuiRunning = 0;
end
%Enter main loop
while count<=OPT.miter,
 %----------------------------------------------------------------------
 % Get cost, prediction error, gradient and Hessian(approx) and save
 % cost into log
 %----------------------------------------------------------------------
 if OPT.cmpgrd,
  [cost,pe]       = feval(FUN,Z,theta,OPT,ARGS,0); 
  cost_log(count) = cost;
  map = eye(length(theta_struct));
 else
  try,
   [cost,pe,grad,phi,map] = feval(FUN,Z,theta,OPT,ARGS,1);
  catch
   [cost,pe,grad,phi]     = feval(FUN,Z,theta,OPT,ARGS,1);
   map = eye(length(theta_struct));
  end
  cost_log(count)        = cost;
 end
 
 %----------------------------------------------------------------------
 % Make sure cost is sensible
 %----------------------------------------------------------------------
 if isnan(cost) | isinf(cost),
  udisp('Cost is either Inf or NaN. Terminating',gui,guih); 
  break; 
 end

 %----------------------------------------------------------------------
 % This part of the code will compute numerical gradients and Jacobians
 % If desired (via OPT.ngt = 1), it will compare the numerical gradient 
 % with the one returned from the user defined function.
 %----------------------------------------------------------------------
 if OPT.ngt || OPT.cmpgrd || OPT.nht,
  epsdiff=(eps)^(2/3);
  npar = length(theta_struct);
  for ngti=1:npar,
   thngt = theta; thngt(theta_struct(ngti)) = thngt(theta_struct(ngti)) + epsdiff;
   if OPT.nht,
    [cngtp,pepos,gradp,phip] = feval(FUN,Z,thngt,OPT,ARGS,1);
   else
    [cngtp,pepos] = feval(FUN,Z,thngt,OPT,ARGS,0);
   end
   thngt = theta; thngt(theta_struct(ngti)) = thngt(theta_struct(ngti)) - epsdiff;
   if OPT.nht,
    [cngtm,peneg,gradn,phin] = feval(FUN,Z,thngt,OPT,ARGS,1);
   else
    [cngtm,peneg] = feval(FUN,Z,thngt,OPT,ARGS,0);
   end
   ngrad(ngti)   = (cngtp-cngtm)/(2*epsdiff);
   nphi(:,ngti)  = (pepos(:)-peneg(:))/(2*epsdiff);
   if OPT.nht,
    nhess(:,ngti) = (gradp(:)-gradn(:))/(2*epsdiff);
   end
  end
  
  %Make sure we copy the lower triangle to the upper
  if OPT.nht,
   for jj=1:npar,
    for ii=jj+1:npar,
     nhess(jj,ii) = nhess(ii,jj);
    end
   end
  end
  
  %If no gradient returned then use numerical gradient
  if OPT.cmpgrd,
   grad = ngrad(:);
   phi  = nphi;
  end
  
  %If user wants to know how close numerical gradient is compared
  %with the one they compute in their own function...
  if OPT.ngt,
   udisp(' ',gui,guih)
   udisp(' ',gui,guih)
   udisp('  Num. Grad.,   Anal. Grad.',gui,guih);
   udisp([ngrad(:) map*grad(:)],gui,guih)
   udisp(sprintf('Numerical Gradient Difference = %13.3e',norm((map*grad(:)-ngrad(:))./(eps+map*grad(:)),'inf')),gui,guih)
   keyboard
  end
 end

 %----------------------------------------------------------------------
 % COMPUTE SEARCH DIRECTION:
 %
 % Try and compute search direction, otherwise just compute neg gradient
 %----------------------------------------------------------------------
 redo = 1;
 while redo
  
  redo = 0;
  switch OPT.dir,
   
   %------------------------------------------------------
   % Robust Gauss-Newton method
   %------------------------------------------------------
   case 'rgn'
    %Get search direction
    [U,S,V] = svd(full(phi),0);
    s       = diag(S);
    loop_count = 0;
    while loop_count<length(s),
     r  = sum(s>=s(1)*svtol);
     s(r+1:end) = s(r+1:end) + svtol./(s(r+1:end)+eps);
     g  = -(V*((U'*pe(:))./s));
      if abs(g'*grad)>=1e-8*max(norm(g)*norm(grad)), break; end
     if svtol<=sqrt(eps),
      OPT.dir = 'trust';
      redo    = 1; %recalculate search direction using new method
     else,
      svtol=max(eps,svtol/2);
     end
     loop_count=loop_count+1;
    end

    %Compute robust norm of gradient
    tnorm = norm(grad,'inf')/max(sqrt(eps),sqrt(cost));
    direction='rGN';
    if strcmp(OPT.dir,'trust'), direction='trust'; end
    
    
    %------------------------------------------------------
    % Gauss-Newton method
    %------------------------------------------------------
   case 'gn'
    %Get search direction
    g = -phi\pe(:); r=size(phi,2);

    %Compute robust norm of gradient
    tnorm = norm(grad,'inf')/max(sqrt(eps),sqrt(abs(cost)));
    direction='GN';

    %------------------------------------------------------
    % Levenberg-Marquardt method
    %------------------------------------------------------
   case 'lm'
    %Get search direction
    [U,S,V] = svd(phi,0); 
    s       = diag(S);
    snew    = s+delta./(s+eps); 
    r       = sum(snew > sqrt(eps));
    g       = -(V(:,1:r)*((U(:,1:r)'*pe(:))./snew(1:r)));
    tnorm = norm(grad,'inf')/max(sqrt(eps),sqrt(abs(cost)));
    direction = 'LM';

    %------------------------------------------------------
    % Trust method (exact solution)
    %------------------------------------------------------
   case 'trust'
    %Hack to invert delta value between LM and trust
    if trust_1st, 
     delta     = 1/OPT.delta;  
     trust_1st = 0;
    end
    [U,S,V] = svd(full(phi),0); 
    s       = diag(S); 
    scal    = 2/(Z.Ny-OPT.n);
    lam     = eps;
    
    %Get trust region delta
    for bisec=1:OPT.lmax,
     sn = s+lam./(s+eps); 
     r  = sum(sn>sn(1)*eps);
     q  = U(:,1:r)'*pe(:);
     p  = -(V(:,1:r)*(q./sn(1:r)));
     if norm(p)>delta,
      hitdel=1;
      for inloop=1:50,
       n1   = p'*p;
       n2   = q'*((q.*s(1:r).^2)./((s(1:r).^2+lam).^3));
       dlam = (n1/n2)*((sqrt(n1)-delta)/delta);
       lam  = max(eps,lam+dlam);
       sn   = s+lam./(s+eps); 
       r    = sum(sn>sn(1)*eps);
       q    = U(:,1:r)'*pe(:);
       p    = -(V(:,1:r)*(q./sn(1:r)));
       if abs(norm(p)-delta)<1e-2*delta, break; end
      end
     else
      hitdel=0;
     end
     
     %tnorm=abs(p'*grad)/max(sqrt(eps),sqrt(abs(cost)));
     tnorm = norm(grad,'inf')/max(sqrt(eps),sqrt(abs(cost)));

     theta_new = theta;
     theta_new(theta_struct)=theta(theta_struct)+map*p; 
     nc=feval(FUN,Z,theta_new,OPT,ARGS,0);
     if isnan(nc), nc=inf; end

     %Compute actual/predicted cost reduction
     
     pp   = phi*p; tmp1=pp'*pe(:)*scal; tmp2=pp'*pp*scal;
     rho  = (cost - nc)/(1*eps-tmp1-0.5*tmp2);

     %Adaptively change trust region according to local performance.
     method = 'new'; %'old' is other option
     if strcmp(method,'new'),
      alp1 = 0.5;  alp2 = 2;    alp3 = 1.01;
      eta1 = 0.01; eta2 = 0.95; eta3 = 1.05;
      if rho <= 0,
       delta = max(alp1*delta,eps);
      elseif rho < eta2,
       delta = delta*(alp1 + (1-alp1)*(rho/eta2)^2);
      elseif rho >= eta2,
       delta = delta*min(1e12,alp3 + (alp2-alp3)*exp(-((rho-1)/(eta2-1))^2));
      end
     else %Old method
      if rho < 0.25,
       delta = max(0.25*norm(p),eps);
      elseif rho > 0.75 & hitdel,
       delta = min(2*delta,1e16);
      end
     end
     
     %If we reduced the cost sufficiently then stop.
     if rho > 0.01, theta=theta_new; break; end
    end
    direction='trust';

    %Store some stuff
    phi_old=phi; grad_old=grad;

    %------------------------------------------------------
    % Trust method (exact solution)
    %------------------------------------------------------
   case 'trust_scaled'
    %Hack to invert delta value between LM and trust
    if count==1, delta=1/delta; end
    lam  = eps; 
    scal = 2/(Z.Ny-OPT.n);
    R    = triu(qr([phi pe]));
    nmp  = size(phi,2);    
    phi  = sqrt(scal)*R(1:nmp,1:nmp);
    pe   = sqrt(scal)*R(1:nmp,end);
    dpi  = diag(phi'*phi);
    r    = nmp;
    
    %Compute a diagonal scaling term.
    %Dscal = diag(min(1/sqrt(eps),max(sqrt(eps),1./(abs(dpi)))));
    
    %Get trust region delta
    for bisec=1:OPT.lmax,
     %Factor phi'*phi + lambda*Dscal
     R = triu(qr([phi;sqrt(lam)*eye(nmp)]));
     R = R(1:nmp,:);
     
     %Solve for p,   R'*Rp=-g
     p = -R\((R')\(grad));
     q = (R')\p;
     
     if norm(p)>delta,
      hitdel=1;
      for inloop=1:50,
       np   = norm(p);
       nq   = norm(q);
       dlam = (np/nq)^2*((np-delta)/delta);
       lam  = max(eps,lam+dlam);

       %Factor phi'*phi + lambda*I
       R = triu(qr([phi;sqrt(lam)*eye(nmp)]));
       R = R(1:nmp,:);
       
       %Solve for p,   R'*Rp=-g
       p = -R\((R')\(grad));
       q = (R')\p;
       
       if abs(norm(p)-delta)<1e-2*delta, break; end
      end
     else
      hitdel=0;
     end
     
     tnorm = norm(grad,'inf')/max(sqrt(eps),sqrt(abs(cost)));

     theta_new = theta;
     theta_new(theta_struct)=theta(theta_struct)+map*p; 
     nc = feval(FUN,Z,theta_new,OPT,ARGS,0);
     if isnan(nc), nc=inf; end

     %Compute actual/predicted cost reduction
     pp   = phi*p; 
     tmp1 = pp'*pe(:); 
     tmp2 = pp'*pp;
     rho1 = (cost - nc);
     rho2 = (1*eps-tmp1-0.5*tmp2);
     
     %Adaptively change trust region according to local performance.
     if rho1 < 0.25*rho2,
      delta = max(0.25*norm(p),eps);
     elseif rho1 > 0.75*rho2 & hitdel,
      delta = min(2*delta,1e16);
     end

     %If we reduced the cost sufficiently then stop.
     if rho1 > 0.01*rho2, theta=theta_new; break; end
    end
    direction='tst_scl';

    %Store some stuff
    phi_old=phi; grad_old=grad;
    
    %------------------------------------------------------
    % BFGS quasi-Newton method
    %------------------------------------------------------
   case 'bfgs'
    %Do BFGS update
    if count>1
     td=theta(theta_struct)-theta_old(theta_struct);
     gd=grad-grad_old;     
     if td'*gd > sqrt(eps)*norm(td)*norm(gd),
      rho=1/(eps+gd'*td);

      %Set initial Hessian inverse for BFGS
      if count==2 && ~hess_supplied, 
       H = Ith*(td'*gd)/(gd'*gd); 
      end
      
      H = (Ith-rho*td*gd')*H*(Ith-rho*gd*td') + rho*td*td';
      
      %me=min(eig(H)); if me<=0, H=H+(2*abs(me)+1e-9)*Ith; end
     end
    elseif ~hess_supplied;
     H = 1/norm(grad)*H;
    end
    
    %Store inverse Hessian approximation
    ARGS.invHess=H;

    %Get search direction
    g=-H*grad;
    tnorm = norm(grad,'inf')/max(sqrt(eps),sqrt(abs(cost)));
    r=length(theta_struct);
    direction='QNLS';

    %------------------------------------------------------
    % BFGS + Trust region quasi-Newton method
    %------------------------------------------------------
   case 'bfgs_trust'
    %Do BFGS update
    if count>1,
     td=theta(theta_struct)-theta_old(theta_struct);
     gd=grad-grad_old;

     if td'*gd > sqrt(eps)*norm(td)*norm(gd),
     
      %Set initial Hessian inverse for BFGS
      if count==2, H = Ith*(gd'*gd)/(td'*gd); end
      
      Htd = H*td;
      H   = H + (gd*gd')/(td'*gd) - (Htd*Htd')/(td'*Htd);
     end
    else
     delta = 1/delta;
     H     = 1/norm(grad)*H;
    end

    ARGS.Hess = H;

    %Use trust region method.
    S=eig(H); 
    if ~isreal(S), 
     udisp('Hessian has complex eigenvalues',gui,guih); 
    end
    lam_min=min(0,min(real(S)));
    r=length(theta_struct);

    if lam_min < 0, lam=-2*lam_min; else, lam_min=0; lam=0; end
    eta=0.01;

    for bisec=0:OPT.lmax,
     %Get search direction
     [cH,pp]=chol(H+lam*eye(size(H)));
     while pp~=0,
      lam=2*(lam+1e-3);
      [cH,pp]=chol(H+(lam+eps)*eye(size(H)));
     end
     p=-(cH\(cH'\grad));

     %Check that search direction satisfies ||dp|| <= del,
     %if not then find Lagrange multiplier for constrained
     %problem.
     if norm(p)>delta,
      hitdel=1;
      for i=1:100,
       q=cH'\p;
       dlam=(norm(p)/norm(q))^2*((norm(p)-delta)/delta);
       lam=max(-lam_min,lam+dlam);
       [cH,pp]=chol(H+(lam+eps)*eye(size(H)));
       while pp~=0,
        lam=2*(lam+1e-1);
        [cH,pp]=chol(H+(lam+eps)*eye(size(H)));
       end
       p=-cH\(cH'\grad);
       if abs(norm(p)-delta) < 1e-3, break; end
      end
     else
      hitdel=0;
     end

     %Get norm of step
     tnorm = norm(grad,'inf')/max(sqrt(eps),sqrt(abs(cost)));

     %Get new cost
     theta_new=theta; 
     theta_new(theta_struct)=theta(theta_struct)+map*p; 
     nc=feval(FUN,Z,theta_new,OPT,ARGS,0);
     if isnan(nc), nc=inf; end

     %Compute actual/predicted cost reduction
     rho = (cost - nc)/(eps-grad'*p-0.5*(p'*(H*p + lam*p)));

     %Adaptively change trust region according to local
     %performance.
     if rho < 0.25,
      delta=0.25*norm(p);
     elseif rho > 0.75 & hitdel,
      delta=min(2*delta,1e16);
     end

     %If delta is too small then quit searching
     if delta<= eps, break; end

     %If we reduced the cost sufficiently then stop.
     if rho > eta, theta_old=theta; theta=theta_new; break; end
    end
    direction='QNTR';

    %------------------------------------------------------
    % Steepest-descent method
    %------------------------------------------------------
   case 'grad'
    %Get search direction
    g=-grad; tnorm=norm(g); if tnorm>eps, g=g/tnorm; end; r=length(theta_struct);
    direction='Grad.';


   otherwise,
    error('OPT.dir not known!');

  end
 end
 %----------------------------------------------------------------------
 % END OF SEARCH DIRECTION COMPUTATION:
 %----------------------------------------------------------------------


 %-----------------------------------------------------------------------
 % COMPUTE STEP LENGTH:
 %
 % Perform bisection determination of step length.
 %----------------------------------------------------------------------
 alp=1; grad_old=grad; panic=0; theta_old = theta;
 c1=0*sqrt(eps); pdn=norm(grad); pdir=grad; if pdn>eps, pdir=pdir/pdn; end
 if ~strcmp(OPT.dir,{'agn','trust','trust_scaled','bfgs_trust'})
  bisec=0;
  while 1,
   if ~panic,
    switch OPT.dir,
     case 'lm',
      sn = s+delta./(s+eps); 
      r  = sum(sn>sqrt(eps));
      p  = -(V(:,1:r)*((U(:,1:r)'*pe(:))./sn(1:r)));      
     otherwise,
      p  = g*alp;
    end
   else
    p=alp*pdir;
   end
   theta_new = theta;
   theta_new(theta_struct) = theta(theta_struct)+map*p;
   nc=feval(FUN,Z,theta_new,OPT,ARGS,0);
   if isnan(nc), nc=inf; end
   if nc<cost+c1*grad'*p, break; end
   alp=alp/2; delta=min(1e16,2*delta); bisec=bisec+1;
   if bisec>=OPT.lmax,
    if panic,
     break;
    else
     direction='Max-BS';
     break;
     panic=1;
     bisec=0;
    end
   end
  end
  if nc<cost, theta=theta_new; else, nc=cost; end
 end
 %-----------------------------------------------------------------------
 % END OF STEP LENGTH CALCULATION VIA BACKTRACKING
 %----------------------------------------------------------------------

 %----------------------------------------------------------------------
 % COMPUTE ADAPTIVE CONDITION NUMBER SVTOL
 %----------------------------------------------------------------------
 if OPT.adapt,
  switch OPT.dir,
   case 'rgn'
    if bisec>5,
     svtol=2*svtol;
    elseif bisec==0,
     svtol=max(sqrt(eps),svtol/4);
    end

   case 'lm'
    if bisec==0, delta=delta/2; end
  end
 end


 %----------------------------------------------------------------------
 % SAVE THINGS
 %----------------------------------------------------------------------
 if OPT.saveit, theta_saved=[theta_saved theta]; end

 %----------------------------------------------------------------------
 % Update ARGS for certain parametrizations (e.g. DDLC)
 %----------------------------------------------------------------------
 if nargin>5, [ARGS,theta]=feval(fargs,theta,ARGS); end


 %----------------------------------------------------------------------
 % Display stuff if meant to
 %----------------------------------------------------------------------
 if OPT.dsp,
  if count==1,
   switch OPT.dir,
    case 'rgn',
     tol_name='SV Tol';
     tol_value=svtol;
    case 'agn',
     tol_name='Delta';
     tol_value=delta;
    case 'gn',
     tol_name='SV Tol';
     tol_value=OPT.subtol;
    case 'lm',
     tol_name='Delta';
     tol_value=delta;
    case 'trust',
     tol_name='Delta';
     tol_value=delta;
    case 'trust_scaled',
     tol_name='Delta';
     tol_value=delta;
    case 'bfgs_trust',
     tol_name='Delta';
     tol_value=delta;
    otherwise
     tol_name='N/A';
     tol_value=0;
   end
   udisp(sglines,gui,guih);
   str1=sprintf('%s%13s%13s%9s%7s%2i)%13s%7s','Iter#','Cost','G-N Norm','Bisec#','SV#(/',size(map,2),tol_name,'Dir');
   udisp(str1,gui,guih);
   udisp(sglines,gui,guih);
   str2=sprintf('%5i%13.3e%13s%9s%10s%13.3e%8s',0,cost,'-','-','-',tol_value,'-');
   udisp([str2],gui,guih)
  end
  switch OPT.dir,
   case 'rgn',
    tol_name='SV Tol';
    tol_value=svtol;
   case 'agn',
    tol_name='Delta';
    tol_value=delta;
   case 'gn',
    tol_name='SV Tol';
    tol_value=OPT.subtol;
   case 'lm',
    tol_name='Delta';
    tol_value=delta;
   case 'trust',
    tol_name='Delta';
    tol_value=delta;
   case 'trust_scaled',
    tol_name='Delta';
    tol_value=delta;
   case 'bfgs_trust',
    tol_name='Delta';
    tol_value=delta;
   otherwise
    tol_name='N/A';
    tol_value=0;
  end
  str1=sprintf('%5i%13.3e%13.3e%9i%10i%13.3e%8s',count,nc,tnorm,bisec,r,tol_value,direction);
  udisp(str1,gui,guih)
 end

 %----------------------------------------------------------------------
 % Update main loop counter
 %----------------------------------------------------------------------
 count=count+1;

 %----------------------------------------------------------------------
 % Check termination conditions (do not test until tried to decrease)
 %----------------------------------------------------------------------
 if count>1,
  %Poll the Stop Estimation Gui Button Pressed variable
  if stopEstGuiRunning,
   stopEst = stopEstGui.buttonPressed;
  else
   stopEst = 0;
  end
  if ~isempty(guih),
   if(guih.gui.eto.myEstModelPopup.check == 1)
    stopEst = 1;
   end
  end
  if stopEst
   %Stopping because the user wants to
   udisp(sglines,gui,guih);
   udisp('Termination due to Users wishes    ',gui,guih);
   udisp(sglines,gui,guih);
   ARGS.whystop='Termination due to Users wishes';
   break;
  elseif isfield(OPT,'sv') & nc < OPT.sv,
   %Termination due to cost falling below given value
   if OPT.dsp,
    udisp(sglines,gui,guih);
    udisp('Termination due to cost falling below given value in OPT.sv    ',gui,guih);
    udisp(sglines,gui,guih);
   end
   ARGS.whystop='Termination due to cost falling below given value in OPT.sv';
   break;
  elseif tnorm<OPT.tol,
   %Termination due to gradient norm less than OPT.tol
   if OPT.dsp,
    udisp(sglines,gui,guih);
    udisp('Termination due to gradient norm less than OPT.tol           ',gui,guih);
    udisp(sglines,gui,guih);
   end
   ARGS.whystop='Termination due to gradient norm less than OPT.tol';
   break;
  elseif count>OPT.miter,
   %Termination due to number of iterations exceeding OPT.miter
   if OPT.dsp,
    udisp(sglines,gui,guih);
    udisp('Termination due to number of iterations exceeding OPT.miter  ',gui,guih);
    udisp(sglines,gui,guih);
   end
   ARGS.whystop='Termination due to number of iterations exceeding OPT.miter';
   break;
  elseif bisec>=OPT.lmax,
   %Termination due to number of bisections exceeding OPT.lmax
   if OPT.dsp,
    udisp(sglines,gui,guih);
    udisp('Termination due to number of bisections exceeding OPT.lmax   ',gui,guih);
    udisp(sglines,gui,guih);
   end
   ARGS.whystop='Termination due to number of bisections exceeding OPT.lmax';
   break;
  elseif delta<2*eps,
   %Termination due to trust region too small
   if OPT.dsp,
    udisp(sglines,gui,guih);
    udisp('Termination due to trust region being too small (< 2*eps)    ',gui,guih);
    udisp(sglines,gui,guih);
   end
   ARGS.whystop='Termination due to number of bisections exceeding OPT.lmax';
   break;
  elseif 0 & abs(cost-nc)<OPT.mdec*abs(p'*grad),
   %Termination due to cost difference less than OPT.mdec
   if OPT.dsp,
    udisp(sglines,gui,guih);
    udisp('Termination due to cost difference less than OPT.mdec        ',gui,guih);
    udisp(sglines,gui,guih);
   end
   break;
  elseif 0 & (bisec<OPT.lmax) & (norm(theta-theta_old)<sqrt(eps)*norm(theta_old)),
   %Termination due to insufficient progress
   if OPT.dsp,
    udisp(sglines,gui,guih);
    udisp('Termination due to parameters not changing                   ',gui,guih);
    udisp(sglines,gui,guih);
   end
   break;
  elseif isinf(delta),
   %Termination due to delta being too big
   if OPT.dsp,
    udisp(sglines,gui,guih);
    udisp('Termination due to delta being too large                     ',gui,guih);
    udisp(sglines,gui,guih);
   end
   break;
  end
 end
end

%If stopEstimation GUI is running then close it
if stopEstGuiRunning,
 stopEstGui.setVisible(false);
 clear stopEstGui;
end

% Make sure we save some important info.
ARGS.OPT   = OPT;
ARGS.theta = theta;
if OPT.miter>0,
 cost_log(end+1) = nc;
 if OPT.saveit, ARGS.thetait = theta_saved; end
else
 nc=feval(FUN,Z,theta,OPT,ARGS,0);
 cost_log(end+1) = nc;
end

if exist('phi','var'),
 ARGS.jacobian = phi;
end

