% STARTH - function to initialise estimate of noise model in case user
% did not specify it.  This function is not meant to be called by a user -
% instead it is just an auxiliary function that is used internally by other
% routines; most importantly EST.m
%
% Usage is:
%
% M = startH(Z,M,OPT)
%
% written by Brett Ninness, School of EE & CS
%                           University of Newcastle
%             		        Australia.

% Copyright (C) Brett Ninness.

function M = startH(Z,M,OPT)

% Extract input and output from data matrix
[y,u,ny,nu,Ny] = Z2data(Z);

%Detect if gui is running
gui = 0; guih = [];
if isfield(OPT,'gui'),
 if ~isempty(OPT.gui)
  gui  = 1;         %GUI is running
  guih = OPT.gui;   %GUI handle
 end
end

if (ny>1 || strcmpi(M.type,'ss'))  % Multiple output | ss model requires SID
 if ~isfield(M,'ss') % Do not overwrite any initial state space system guess
  if OPT.dsp 
   udisp('Finding initial noise Model by SID...',gui,guih); 
  end;
  M = sid(Z,M);
  
 else % OK, some ss bits are specified, but all of them?
  if ~isfield(M.ss,'A'),
   error('Bits of M.ss pre-specified, but not enough - M.A..?');
  else,
   nx=size(M.ss.A,1);
  end
  if ~isfield(M.ss,'R') M.ss.R=0.01*eye(ny,ny);      end;
  if ~isfield(M.ss,'Q') M.ss.Q=100*eye(nx,nx);      end;
  if ~isfield(M.ss,'C') M.ss.C=randn(ny,nx);    end;
 end
else  % Scalar output case
 switch M.type,
  
  case {'oe','fir'}
   M.C=1; M.D=1;  % In that case, noise model is fixed at H=1;
   
  case {'arx'}
   M.C=1; M.D=1;  % Noise model is H=1/A, which comes from Ay = Bu + e - NOT D;
   
  case {'arma'}
   if OPT.dsp
    udisp('Finding initial noise Model by Hannan-Rissanen...',gui,guih); 
   end
   % Fit high order AR model to y
   MM       = M; 
   MM.A     = min(100,floor(Ny/10)); 
   MM.nA    = MM.A; 
   MM.type  ='ar'; 
   MM.op    = 'q';
   MM.delay = 0;
   opt      = OPT; 
   opt.fast = 1;%  Make initialisation fast
   g        = barx(y,MM,opt);
   
   % So now we have model y = 1/Ae => estimate e as Ay:
   e = filter(g.A,1,y);  
   e = e(length(g.A):length(e)); 
   v = y(length(g.A):length(y));
   
   % Finally, estimate C & D in ARX model Dv = Ce:
   MM       = M;
   MM.A     = MM.nA;
   MM.B     = MM.nC;
   MM.nB    = MM.nC;
   MM.nu    = 1;
   MM.type  = 'arx'; 
   MM.delay = 0;
   
   g   = barx([v(:),e(:)],MM,opt);
   M.A = g.A;
   M.C = g.B; 
   M.C = M.C/M.C(1);  
   % Check model makes sense in that it is stable and minimum phase
   if M.op=='q'
    M.C = stab(M.C,'q'); 
    M.A = stab(M.A,'q');
   else
    M.C = stab(M.C,'d',M.T); 
    M.A = stab(M.A,'d',M.T);
   end; 
   
  case {'armax'}
   if OPT.dsp
    udisp('Finding initial noise Model by Hannan-Rissanen...',gui,guih);
   end
   % Get estimate of signal (B/A)u via current guess at B/A.
   nu = M.nu;
   w  = zeros(size(y));
   if (M.op=='q')
    for k=1:nu
     len = max(length(M.A),length(M.B(k,:)))-1; % First len samples of w forced to match those of y;
     w   = w+filter(M.B(k,:),M.A,u(:,k));%,filter(M.A(k,:),1,y(1:len))-filter(M.B(k,:),1,u(1:len,k)));
    end;
   else
    for k=1:nu
     w = w+delfilter(M.B(k,:),M.A,u(:,k),M.T);%,zeros(1,length(M.A(k,:))-1));
    end;
   end;

   % This allows estimate of coloured noise v = y - Gu: Fit high order AR model to it.
   MM       = M; 
   MM.A     = min(100,floor(Ny/10)); 
   MM.nA    = MM.A; 
   MM.B     = 1;  
   MM.nB    = 1; 
   MM.type  = 'ar'; 
   MM.op    = 'q';
   MM.delay = 0;
   MM.nu    = 0;
   opt      = OPT; 
   opt.fast = 1; 
   v        = y(:)-w(:); 
   g        = barx(v,MM,opt); 
   
   % So now we have model w = 1/De => estimate e as Dw:
   e = filter(g.A,1,v);  
   e = e(length(g.A):length(e)); 
   v = v(length(g.A):length(v));
   
   % Finally, estimate C & D in ARX model Dv = Ce:
   MM       = M; 
   MM.A     = []; 
   MM.B     = M.C; 
   MM.nA    = 0; 
   MM.nB    = M.nC;
   MM.op    = M.op;  
   MM.type  = 'arx'; 
   MM.delay = 0;
   MM.nu    = 1;
   v = filter(1,M.A(1,:),v);  % In armax case denominator of noise model fixed as D=A;
   g = barx([v(:),e(:)],MM,opt);
   M.C = g.B; 
   M.C = M.C/M.C(1); 
   
   % Check model makes sense in that it is stable and minimum phase
   if M.op=='q'
    M.C = stab(M.C,'q'); 
   else
    M.C = stab(M.C,'d',M.T); 
   end
   
  case {'bj'}
   if (M.nC>M.nD) M.nC=max(M.nC,M.nD); M.nD=max(M.nC,M.nD);  end;  % Bug catching
   if OPT.dsp 
    udisp('Finding initial noise Model by Hannan-Rissanen...',gui,guih); 
   end
   % Get estimate of signal (B/A)u via current guess at B/A.
   w  = zeros(size(y));
   nu = M.nu;
   np = M.nA;
   if (M.op=='q')
    for k=1:nu
     %len = max(length(M.A(k,:)),length(M.B(k,:)))-1; % First len samples of w forced to match those of y;
     w = w+filter(M.B(k,:),M.A(k,:),u(:,k));%,filter(M.A(k,:),1,y(1:len))-filter(M.B(k,:),1,u(1:len,k)));
    end;
   else
    for k=1:nu
     w = w+delfilter(M.B(k,:),M.A(k,:),u(:,k),M.T);%,zeros(1,length(M.A(k,:))-1));
    end;
   end;
   
   % This allows estimate of coloured noise v = y - Gu: Fit high order AR model to it.
   MM       = M; 
   MM.A     = min(100,floor(Ny/10)); 
   MM.nA    = MM.A; 
   MM.B     = 1;  
   MM.nB    = 1; 
   MM.type  = 'ar'; 
   MM.op    = 'q';
   MM.nu    = 0;
   opt      = OPT; 
   opt.fast = 1; 
   v        = y(:)-w(:); 
   g        = barx(v,MM,opt); 
   
   % So now we have model w = 1/De => estimate e as Dw:
   e = filter(g.A,1,v);  e=e(length(g.A):length(e)); v=v(length(g.A):length(v));
   % Finally, estimate C & D in ARX model Dv = Ce:
   MM       = M; 
   MM.A     = M.nD; 
   MM.B     = M.nC; 
   MM.nA    = M.nD; 
   MM.nB    = M.nC;
   MM.op    = M.op;  
   MM.type  = 'arx'; 
   MM.delay = 0;
   MM.nu    = 1;
   g        = barx([v(:),e(:)],MM,opt);
   M.C      = g.B; 
   M.C      = M.C/M.C(1);
   M.D      = g.A; 
   % Check model makes sense in that it is stable and minimum phase
   if M.op=='q'
    M.C = stab(M.C,'q'); 
    M.D = stab(M.D,'q');
   else
    M.C = stab(M.C,'d',M.T); 
    M.D = stab(M.D,'d',M.T);
   end
 end
end
