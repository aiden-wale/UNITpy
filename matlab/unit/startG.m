% STARTG - function to initialise estimate of dynamics model in case user
% did not specify it.  This function is not meant to be called by a user -
% instead it is just an auxiliary function that is used internally by other
% routines; most importantly EST.m
%
% Usage is:
%
% M = startG(Z,M,OPT)
%
% written by Brett Ninness, School of EE & CS
%                           University of Newcastle
%                       Australia.

% Copyright (C) Brett Ninness.

function M = startG(Z,M,OPT)

% Extract input and output from data matrix
OPT = startOPT(OPT);
Z   = startZ(Z);
[y,u,ny,nu,Ny] = Z2data(Z);

% Check that model supports data 
[flag,message] = checkZM(Z,M);
if flag, error(message); end

% Make sure input, output and state dimension come from the model.
nu = M.nu;
ny = M.ny;
nx = M.nx;

% Detect if gui is running
gui = 0; guih = [];
if isfield(OPT,'gui'),
 if ~isempty(OPT.gui)
  gui  = 1;         %GUI is running
  guih = OPT.gui;   %GUI handle
 end
end

if nu>0 % Only initialise I/O dynamics if there are some to estimate!

 % OK, now initialise estimate of dynamics if an initial estimate has not been specified
 % Firstly, if multi-output, ss or bilinear model, then must use subspace method
 if (M.ny>1 | any(strcmpi(M.type,{'ss','bilin','bilinear'})))
  if ~isfield(M,'ss') % Do not overwrite any initial state space system guess
   opt=OPT; opt.alg='n4sid';
   switch Z.type
    case 'time'
     if OPT.dsp,
      udisp('Finding initialisation for estimate of Dynamics via Subspace ID...',gui,guih)
     end
     if Z.T==0 || M.op=='s',
      M.T    = mean(diff(Z.t));
      M.op   = 'q';
      m      = sid(Z,M,OPT);
      sys    = ss(m.ss.A,m.ss.B,m.ss.C,m.ss.D,m.T);
      sys    = d2c(sys,'tustin');
      m.ss.A = sys.A; 
      m.ss.B = sys.B; 
      m.ss.C = sys.C;
      if M.estK,
       try
        [junk1,junk2,K]=care(m.ss.A',m.ss.C',eye(M.nx),1000*eye(M.ny));
        m.ss.K = K';
       catch
        m.ss.K = zeros(M.nx,M.ny);
       end
      else
       m.ss.K = [];
      end
      m.estD = 0; %we do not allow for a D term in continuous model
      m.ss.D = []; 
      if M.estX1,
       m.ss.X1 = zeros(M.nx,1);
      else
       m.ss.X1 = [];
      end
      m.op   = 's'; 
      m.T    = Z.T;
     else
      m = subspace(Z,M,OPT);
     end
    case 'frequency'
     if OPT.dsp,
      udisp('Finding initialisation for estimate of Dynamics via Freq. Dom. Subspace ID...',gui,guih)
     end
     m = subspace(Z,M,OPT);
    otherwise
     error('Data type (Z.type) not known');
   end
   m.type=M.type; M=m; % Make sure (possible bilinear type) is not overwritten
   if any(strcmpi(M.type,{'bilin','bilinear'})),
    M.ss.F = zeros(m.nx,m.nx*m.nu); M.ss.G = zeros(m.ny,m.nx*m.nu);
   end
  else % OK, some ss bits are specified, but all of them?
   nx = M.nx; nu = M.nu; ny = M.ny;
   if isfield(M,'ss'),
    if ~isfield(M.ss,'A'),  M.ss.A  = zeros(nx);       end
    if ~isfield(M.ss,'B'),  M.ss.B  = zeros(nx,nu);    end
    if ~isfield(M.ss,'C'),  M.ss.C  = zeros(ny,nx);    end
    if ~isfield(M.ss,'D'),  M.ss.D  = zeros(ny,nu);    end
    if ~isfield(M.ss,'Q'),  M.ss.Q  = eye(nx)/10;    end
    if ~isfield(M.ss,'S'),  M.ss.S  = zeros(nx,ny);    end
    if ~isfield(M.ss,'R'),  M.ss.R  = eye(ny)/1000;    end
    if ~isfield(M.ss,'X1'), M.ss.X1 = zeros(nx,1);     end
    if any(strcmpi(M.type,{'bilin','bilinear'})),
     if ~isfield(M.ss,'F'),  M.ss.F  = zeros(nx,nx*nu); end
     if ~isfield(M.ss,'G'),  M.ss.G  = zeros(ny,nx*nu); end
    end
   end
   if ~isfield(M.ss,'K')
    % Compute Kalman Gain for innovations form representation
    try,
     [P,dd,K] = dare(M.ss.A',M.ss.C',M.ss.Q,M.ss.R);
     M.ss.K   = K';
    catch
     K        = zeros(nx,ny);
    end
   end;
  end;

 
 else  %We are looking at a transfer function model

  %Is it time or frequency domain data?
  switch Z.type,

   case 'time',
    %For each input, check if an initial parameter
    %guess was supplied and estimate one if not. This is done based on the
    %order information in M.nA and M.nB
    
    %  Make initialisation fast
    opt      = OPT; 
    opt.fast = 1; 
    
    Msv = M;
    
    %Must handle BJ and OE models separately because they allow
    %a different A poly for each input, this means that we must
    %estimate a SISO model for each input. For all other cases,
    %we can estimate just the one A polynomial.
    switch M.type
     case {'bj','oe'}
      %Loop over inputs
      for ku = 1:M.nu,
       %Determine if we should estimate for this pair
       %based on what was supplied in M.A and M.B
       estimate = 0; %Initially don't estimate anything
       if M.nA(ku) ~= size(M.A,2)-1,
        estimate = 1;
       end
       if M.nB(ku) ~= size(M.B,2)-1,
        estimate = 1;
       end
       
       if estimate,
        if OPT.dsp
         udisp('Finding initial dynamics model via Steiglitz-McBride...',gui,guih);
        end
            
        % Get initial ARX estimate of order specified for i/o model we are up to
        ZZ    = Z;
        ZZ.u  = ZZ.u(:,ku);
        ZZ.nu = 1;
        MM    = M;
        MM.nB = MM.nB(ku);
        MM.B  = MM.nB;
        MM.nA = MM.nA(ku);
        MM.A  = MM.nA;
        MM.nu = 1;
        MM.ny = 1;
        MM.op = 'q';  % M.op could be 's', but we definitely want MM.op='q';
        
        % If we are getting an initial estimate for an eventual
        % CT model, we get better results now with relative degree zero
        if (M.op=='s')   
         MM.nB = MM.nA;   
        end;
        
        %Now call barx to get estimate
        g       = barx(ZZ,MM,opt);
        gsave   = g;
           
        %Make sure the A polynomial is stable
        g.A = stab(g.A,g.op,g.T);
        
        %Now save the cost function value so that we
        %can determine which iteration of the SM method
        %gives best results
        mm      = g;
        mm.type = 'oe';
        mm.C    = 1;
        mm.D    = 1;
        costold = VN(ZZ,m2theta(mm),opt,mm,0);
        mi      = 1;
        
        if OPT.dsp,
         udisp(sprintf('it#: %3i   MSE cost: %3.5e',0,costold),gui,guih);
        end
        
        % Then prefilter and estimate again to minimize noise induced bias
        % error and also implicit high frequency distortion in fit.
        for k=1:OPT.smits % Usually only once around is needed, 4 => robustness to rare cases
         if (M.op=='q')
          yf = filter(sum(g.A),g.A,ZZ.y);
          uf = filter(sum(g.A),g.A,ZZ.u);
         else
          yf = delfilter(g.A(end),g.A,ZZ.y,M.T);
          uf = delfilter(g.A(end),g.A,ZZ.u,M.T);
         end;
         % Could have estimated an unstable system in this first pass - check for this
         if any(isnan(yf)) error('Unstable initialisation generated - initialise another way'); end;
         MM.delay = 0;
         g        = barx([yf,uf],MM,opt);
         
         %Make sure the A polynomial is stable
         g.A = stab(g.A,g.op,g.T);
         
         %Get new cost
         mm      = g;
         mm.type = 'oe';
         mm.C    = 1;
         mm.D    = 1;
         costnew = VN(ZZ,m2theta(mm),opt,mm,0);
         
         if OPT.dsp,
          udisp(sprintf('it#: %3i   MSE cost: %3.5e',k,costnew),gui,guih);
         end
         
         if costnew<costold, costold=costnew; gsave=g; mi=k; end
        end
        
        %Select best initial point out of gsave according to lowest cost function
        g = gsave;
        if OPT.dsp
         udisp(sprintf('Best results achieved using %i iteration(s)',mi),gui,guih)
        end
        if ku==1,
         Msv = M;
         mna = max(Msv.nA);
         mnb = max(Msv.nB);
         Msv.A = [g.A zeros(1,mna-Msv.nA(ku))];
         Msv.B = [g.B zeros(1,mnb-Msv.nB(ku))];
        else
         Msv.A = [Msv.A; [g.A zeros(1,mna-Msv.nA(ku))]];
         Msv.B = [Msv.B; [g.B zeros(1,mnb-Msv.nB(ku))]];
        end
       end  % Test on whether we need to estimate or not for initial value
      end   % For loop over number of inputs 
      M = Msv;
      
      if strcmp(M.op,'s') 
       % Heads up - we've got a live one - continuous time TF 
       % model from time domain data wanted - achieve via grey box ss
   
       M.theta   = [M.B(1:end)';M.A(2:end)'];
       M.t2m     = 't2m_soe';
       M.type    = 'ss';
       M.par     = 'grey';
       M.alg     = 'gn';
       M.finishM = 'finishMctstf';
       M = feval(M.t2m,M,M.theta);
       
       if estimate   % If this model was estimated, it was a DT one
         nB = M.nB; M.nB = M.nA;     % I also had zero relative degre
         M = feval(M.t2m,M,M.theta); % Convert to ss form
         M.T = mean(diff(Z.t));      % Could be irregularly sampled
         M = d2c(M); % Convert the DT model to a CT initial one
         M.nB = nB;  % Recall the orginal user specified relative degree
         M.B = M.B(end-M.nB:end);  % CT model rel degree not necessarily 1
         M.theta   = [M.B(1:end)';M.A(2:end)'];
       end; 
      end; 
      
     otherwise % M.type must be ar, arma, arx, armax, fir
      %Determine if we should estimate for this pair
      %based on what was supplied in M.A and M.B
      estimate = 0; %Initially don't estimate anything
      if M.nA ~= size(M.A,2)-1,
       estimate = 1;
      end
      if max(M.nB) ~= size(M.B,2)-1,
       estimate = 1;
      end
       
      if estimate,
       if OPT.dsp
        udisp('Finding initial dynamics model via Steiglitz-McBride...',gui,guih);
       end    
       % Now call barx to get estimate
       g       = barx(Z,M,opt);
       gsave   = g;
            
       %Now save the cost function value so that we
       %can determine which iteration of the SM method
       %gives best results
       mm      = g;
       mm.type = 'oe';
       mm.C    = 1;
       mm.D    = 1;
       mm.A    = mm.A(ones(M.nu,1),:);
       costold = VN(Z,m2theta(mm),opt,mm,0);
       mi      = 1;
       
       % Then prefilter and estimate again to minimize noise induced bias
       % error and also implicit high frequency distortion in fit.
       for k=1:OPT.smits % Usually only once around is needed, 4 => robustness to rare cases
        if (M.op=='q')
         yf = filter(sum(g.A),g.A,Z.y);
         uf = filter(sum(g.A),g.A,Z.u);
        else
         yf = delfilter(g.A(end),g.A,Z.y,M.T);
         for knu=1:M.nu
          uf(:,knu) = delfilter(g.A(end),g.A,u(:,knu),M.T);
         end;
        end;
        % Could have estimated an unstable system in this first pass - check for this
        if any(isnan(yf)) error('Unstable initialisation generated - initialise another way'); end;
        MM       = M;
        MM.delay = 0;
        g        = barx([yf,uf],M,opt);
        mm       = g;
        mm.type  = 'oe';
        mm.C     = 1;
        mm.D     = 1;
        mm.A     = mm.A(ones(M.nu,1),:);
        costnew  = VN(Z,m2theta(mm),opt,mm,0);
        
        if costnew<costold, costold=costnew; gsave=g; mi=k; end
       end
        
       %Select best initial point out of gsave according to lowest cost function
       g = gsave;
       if OPT.dsp
        udisp(sprintf('Best results achieved using %i iteration(s)',mi),gui,guih)
       end
       M.A = g.A;
       M.B = g.B;
      end
    
    end  % end switch over model types

   case 'frequency'
    % Extract out relevant vectors from input data
    [F,w,ny,nu,Ny] = Z2data(Z); F=squeeze(F); F=F(:); M.wmax = max(w);
    
    %  Establish frequency domain variable appropriate to time domain operator
    if (M.op=='q'), 
     ww = exp(j*M.w*M.T);
    elseif (M.op=='d'), 
     ww = (exp(j*M.w*M.T)-ones(size(M.w)))/M.T;
    else
     ww = j*M.w; 
    end

    %  Is frequency normalisation necessary?
    M.normw=0;  
    if [M.op == 's', OPT.basis ~= 'ortho'], 
     M.normw = 1; 
    end

    % Check to see of only integer orders where specified as initial guesses
    % for dynamics: if so get initial estimate by fitting ARX model structure.

    if [length(M.A)<2  floor(M.A(:)')==M.A(:)']
     % Get initial ARX estimate of specified order;
     g     = farx(Z,M,OPT); 
     M.th0 = g.th; 
     M.X   = g.X;
     M.A   = g.A; 
     M.B   = g.B; 
     M.n   = length(g.B);
    else  % Otherwise, re-express initial guess wrt chosen basis
     ff = polyval(M.B,ww)./polyval(M.A,ww);  % Response of initial guess
     % Use farx to translate this initial guess to requested basis.
     g     = farx([ff(:),M.w(:)],M,OPT); 
     M.th0 = g.th; 
     M.X   = g.X;
     M.A   = g.A; 
     M.B   = g.B; 
     M.n   = length(g.B);
    end;
    M.theta = M.th0;

  end %End of switch on data type.
 end; % End of test on multiple output | ss type
else
 M.B = 0.0; 
 M.A = 1.0;
end;









