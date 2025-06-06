% This function calls underlying gradient-based search algorithms, based on
% the type of model and type of data passed in. As new models and data
% types are incorporated into the toolbox, this function would expand to
% cater for these new categories if appropriate.
%
% The procedure is to identify which model type is being used
% (M.type), and then to identify what data type is available (Z.type), and
% then to set the cost_fun variable to point to the corresponding cost
% function for this model/data combination.
%
%    Written by Brett Ninness,  School of EE & CS
%               Adrian Wills    University of Newcastle
%                     		    Australia.

% Copyright (C) Brett Ninness

function G = gn(Z,M,OPT);

% Make doublely sure Z is OK
[y,u,ny,nu,Ny,Z] = Z2data(Z);

% If we have a state space model and we are asked not to estimate something,
% then we should take it out now.
if isfield(M,'estD'),  if ~M.estD,  M.ss.D  = []; end; end
if isfield(M,'estK'),  if ~M.estK,  M.ss.K  = []; end; end
if isfield(M,'estF'),  if ~M.estF,  M.ss.F  = []; end; end
if isfield(M,'estG'),  if ~M.estG,  M.ss.G  = []; end; end
if isfield(M,'estX1'), if ~M.estX1, M.ss.X1 = []; end; end
  
% Take parameters in model structure that cannot be found in closed
% form, and stack them into a parameter vector.
theta = m2theta(M);

% Unspecified parts of regularisation model -> defaults
if isfield(OPT,'M'),
 OPT.M = startM(Z,OPT.M);
else
 OPT.M = theta2m(theta*0,M,1);
end

%Detect if gui is running
gui = 0; guih = [];
if isfield(OPT,'gui'),
 if ~isempty(OPT.gui)
  gui  = 1;         %GUI is running
  guih = OPT.gui;   %GUI handle
 end
end

% Save operator in case we change it due to bilinear transform
M.opsv = M.op;

switch M.type,

 %----------------------------------------------------------------------
 %  SS or BILINEAR
 %----------------------------------------------------------------------
 case {'ss','bilin','bilinear'},
  %Determine cost from Z.type
  switch Z.type,
   case 'time',
    if Z.T==0,
     M.costfcn = 'VNsstv';  %Data is irregularly spaced, needs time-varying predictor
    else
     M.costfcn = 'VNss';    %Data is regularly spaced
    end

   case 'frequency',
    M.costfcn='VNssf';  %Data is frequency domain

    if strcmp(M.op,'s') & strcmp(OPT.smeth,'bilin'),
     Z.wsv  = Z.w;
     T      = 2*pi/max(Z.w);
     Z.w    = 2*atan(Z.w*T/2);
     M.Tsv  = M.T;
     M.T    = 1;
     n      = size(M.ss.A,1);
     M.ss.A = pinv(eye(n)-(T/2)*M.ss.A)*((T/2)*M.ss.A+eye(n));
     M.ss.B = (sqrt(T)/2)*(M.ss.A*M.ss.B+M.ss.B);
     M.ss.C = (sqrt(T)/2)*(M.ss.C*M.ss.A+M.ss.C);
     M.ss.D = M.ss.D+M.ss.C*pinv(eye(n)+M.ss.A)*M.ss.B;
     M.op   = 'q';
     if OPT.dsp,
      udisp('Using Bilinear transform to handle continuous domain data.',gui,guih)
     end
    end

   otherwise,
    error('Z.type not specified');
  end

  % If M.par is struct then we must store the structure in M.theta_struct as follows
  if strcmpi(M.par,'struct'),
   M.theta_struct = find([M.ss.Ai(:);M.ss.Bi(:);M.ss.Ci(:);M.ss.Di(:);M.ss.Ki(:);M.ss.X1i(:)]);
  end
  
  % Find initial cost
  cost0 = feval(M.costfcn,Z,theta,OPT,M,0);
  
  % Can only proceed if initial guess implies stable predictor
  if (~isnan(cost0) & (cost0<inf)),
   
   % Do the search
   [theta,cost,G] = argmin(Z,M.costfcn,theta,OPT,M);

   % Convert from stacked vector form -> model structure form
   G = theta2m(theta,G,0);

   % Reverse bilinear transform if necessary
   if strcmp(Z.type,'frequency') & strcmp(M.opsv,'s') & strcmp(OPT.smeth,'bilin'),
    aa     = inv(eye(size(G.ss.A))+G.ss.A);
    G.ss.A = (2/T) * aa * (G.ss.A-eye(size(aa)));
    G.ss.D = G.ss.D - G.ss.C*aa*G.ss.B;
    G.ss.B = (2/sqrt(T))*aa*G.ss.B;
    G.ss.C = (2/sqrt(T))*G.ss.C*aa;
    G.op   = M.opsv;
    Z.w    = Z.wsv;
    G.T    = M.Tsv;
   end

   G.mse = [cost0,cost];      % Output evolution of mean square cost.
   G.var  = cost(end);        % Terminal cost is estimate of innovations variance

   % Because tf model not directly estimated
   G = sstotf(G);

   G.alg='gn';

  else  % If initial predictor is not stable:
   cost0=[]; cost=inf;
   udisp('--------------------------------------------------',gui,guih);
   udisp('No iteration because starting point was unstable !',gui,guih);
   udisp('--------------------------------------------------',gui,guih);
  end

  %----------------------------------------------------------------------
  %  ARMA, ARMAX, OE, BJ,
  %  (and non-linear versions), NARX, NFIR, NARMA, NARMAX, NOE, NBJ
  %----------------------------------------------------------------------
 
 case {'arma','armax','oe','bj','arx','fir'},
  %Determine cost from Z.type
  switch Z.type,
   case 'time',
    if Z.T==0,
     error('Currently do not cater for polynomial models and continuous data!');
     % NOTE - this could be handled via VNcss if in addition
     % it handled structured state-space models - TODO.
    else
     M.costfcn='VN'; % Data is time discrete
    end

   case 'frequency',
    M.costfcn='VNf';

   otherwise,
    error('Data type (Z.type) not known');
  end

  cost0 = feval(M.costfcn,Z,theta,OPT,M,0);  % Check validity of initial estimate.

  if (~isnan(cost0) & (cost0<1e200)),
   % Do the estimation via damped GN line search
   [theta,cost,M] = argmin(Z,M.costfcn,theta,OPT,M);
   
   % Convert from stacked vector form -> model structure form
   G = theta2m(theta,M);

   G.mse = [cost0,cost];      % Output evolution of mean square cost.

   % Terminal prediction error provides estimate of innovations variance
   [costf,pef] = feval(M.costfcn,Z,theta,OPT,G,0);
   G.var  = pef(:)'*pef(:)/(Z.Ny-OPT.n);

   % Because ss model not directly estimated,
   G = tftoss(G);

   % Switch according to data type,
   switch Z.type
    case 'time'

     % May only have estimated nonlinearity, in this case get linear part now
     if strcmpi(M.type,'arx'),
      M.in  = G.in;
      M.out = G.out;
      G     = barx(Z,M,OPT);
     elseif strcmpi(M.type,'fir'),
      M.in  = G.in;
      M.out = G.out;
      G     = onid(Z,M,OPT);
     else  % If not a closed form job, then need to calculate freq resp etc.
      % Calculate error bounds unless doing doing fast version or using ss mod struc.
      if [~OPT.fast,~strcmpi(M.type,'ss')],
       % Get Estimate of white noise variance by sample variance of residuals
       [cost,pe,grad,R] = feval(M.costfcn,Z,theta,OPT,M,1);
       R=(R'*R);  % Need to square up because R is square-root of Hessian.
       G.P = G.var*pinv(R);  % Use pinv since R may be singular when model over-parameterised
       G.th = theta;
      end;  % Check on OPT.fast

      G.alg='gn'; % Record that Gauss-Newton search algorithm was used
     end; % Check on whether linear part was found iteratively

     % Finally, fir/arx with non-linearity implies cov of non-lin bits have to be appended to cov of lin bits
     if [~OPT.fast (strcmpi(M.type,'arx') | strcmpi(M.type,'fir')) ]
      % Get covariance matrix on non-linear component estimates
      [cost,pe,grad,R] = VN(Z,theta,OPT,M,1);
      R=(R'*R);  %Need to square up because R is square-root of Hessian.
      Peta=G.var*pinv(R);
      % Augment with covariance of linear part;
      [x1,x2]=size(G.P); [x2,x3]=size(Peta); G.P = [G.P,zeros(x1,x3);zeros(x2,x1),Peta];
     end;

    case 'frequency'
     % Pack results into output data structure.
     G.delay = M.delay;
     G.T     = M.T;
     G.w     = M.w;
     G.op    = M.op;
     G.th    = theta;
     G.type  ='oe';
     G.C     = [];
     G.D     = [];

     % Add legend for prospective plotting
     G.disp.legend=['Estimated ',G.type,' model using',G.op,' operator'];

     G.alg='gn'; % Record that Gauss-Newton search was employed
   end
  else
   udisp('--------------------------------------------------',gui,guih);
   udisp('No iteration because starting point was unstable !',gui,guih);
   udisp('--------------------------------------------------',gui,guih);
  end;

  %----------------------------------------------------------------------
  %  STATIC
  %----------------------------------------------------------------------
 case {'static'},
  M.costfcn = 'VN';

  cost0 = feval(M.costfcn,Z,theta,OPT,M,0);  % Check validity of initial estimate.

  if (~isnan(cost0) & (cost0<1e200)),
   % Do the estimation via damped GN line search
   [theta,cost,M] = argmin(Z,M.costfcn,theta,OPT,M);
   
   % Convert from stacked vector form -> model structure form
   G = theta2m(theta,M);

   G.mse = [cost0,cost];      % Output evolution of mean square cost.

   % Terminal prediction error provides estimate of innovations variance
   [costf,pef] = feval(M.costfcn,Z,theta,OPT,G,0);
   G.var  = pef(:)'*pef(:)/(Z.Ny-OPT.n);
  else
   udisp('--------------------------------------------------',gui,guih);
   udisp('No iteration because starting point was unstable !',gui,guih);
   udisp('--------------------------------------------------',gui,guih);
  end
  
 otherwise,
  error('Value in M.type is not known!');
end

% Add legend for prospective plotting
G.disp.legend=['Estimated ',upper(G.type),' model:',G.op,' operator:GN search'];
