%   GERROR.M This routine provides a translation from parameter uncertainties in
%   estimated dynamic models G=B/A to 95% confidence ellipses for the
%   estimated frequency response. Also, 95% confidence regions for estimated
%   magnitude and phase are generated
%
%   This function is not intended to be called directly by a user -
%   rather it is in `internal' routine that is used by estimation
%   functions such as est.m, barx.m and onid.m.
%
%   Usage is:
%
%   [Ge,Gvar] = gerror(G)
%
%   where
%
%   G.A      = Denominator estimate polynomial(s).
%   G.B      = Numerator estimate polynomial(s).
%   G.P      = Covariance Matrix of Estimated Parameters.
%   G.T      = Sampling period in seconds. Default is 1s.
%   G.w      = Vector of frequencies (in rad/s *not* normalised frquency)
%              at which to supply estimated frequency response.
%              You only need to supply this if Ghat is requested.
%   G.op     = set to 'q' for shift and 'd' for delta.
%              Default is 'q' if not specified.
%   Ge       = Matrix defining confidence regions.  In the case of MISO
%              systems, this matrix is three dimensional with each `page'
%              representing the error bounds for one input-output model.
%   Gvar     = vector which is var(G(w)), one element per element in G.w.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.


% Copyright (C) Brett Ninness.

function [Ge,Gvar] = gerror(G)

pn    = 20;             % Number of points on confidence ellipse
level = 6;              % 6=>95% Confidence region
nell  = 50;             % Maximum number of confidence region ellipses to plot.
[nu,dummy] = size(G.A); % nu = number of inputs in model structure.

% Check to see how denominator was normalised.
if (G.A(1,1) == 1) lnorm = 1; else lnorm = 0; end;

% Get appropriate frequency domain variable depending on operator used
if (G.op=='q') ww = exp(-j*G.w*G.T); else ww = (G.T)./(exp(j*G.w*G.T)-ones(size(G.w))); end; ww = ww(:);

for r=1:nu 
	Aw(:,r)   = polyval(fliplr(G.A(r,:)),ww); 
	aG(1,r,:) = unwrap(angle(G.G(1,r,:))); 
end
kinc = ceil(length(G.w)/nell);   % Plot no more than `nell' confidence ellipses.

% Initialise matrices used to store error bound information
Ge = zeros(4+pn,floor(length(G.w)/kinc),nu);
mupper = zeros(1,floor(length(G.w)/kinc));  % Initialisation of arrays
mlower = zeros(1,floor(length(G.w)/kinc));  % for storing results.
pupper = zeros(1,floor(length(G.w)/kinc));
plower = zeros(1,floor(length(G.w)/kinc));

Bind = 1; % Which block in theta pertaining to G.B we are up to
Aind = 0; % Which block in theta pertaining to G.A we are up to
for r=1:nu  %  One set of error bounds for each input-output model
  for k=1:length(G.w)
   % GAMMA is derivative of freq. resp. wrt the parameters
%    if strcmp(lower(G.type),'fir')  % In FIR case GAMMA is already computed by onid
%      index = r:nu:(length(G.poles)+length(find(imag(G.poles)))-1)*nu+r;
%      GAMMA = G.GAMMA(k,:); P = G.P(index,index);
%    else  % Otherwise need to figure out what GAMMA is
    if 1,
    kk = 1:1:G.nA(r); kkk = 0:1:G.nA(r)-1; ll = 0:1:G.nB(r);
    if lnorm dGd = (ww(k).^kk)/Aw(k,r); else GAMMA = (ww(k).^kkk)/Aw(k,r); end;
    dGn = (ww(k).^ll)/Aw(k,r);
    if strcmp(lower(G.type),'arx')  % This is special in MISO case since only one den estimationed
     GAMMA = [dGn,-G.G(1,r,k)*dGd];
     index = Bind:Bind+G.nB(r); index = [index,sum(G.nB)+nu+1:sum(G.nB)+nu+G.nA(1)];
     P=G.P(index,index);
    else
     GAMMA = [zeros(1,Bind-1),dGn,zeros(1,sum(G.nB)+nu-G.nB(r)-Bind),...
     zeros(1,Aind),-G.G(1,r,k)*dGd,zeros(1,sum(G.nA)-G.nA(r)-Aind)];
     P=G.P(1:sum(G.nB)+nu+sum(G.nA),1:sum(G.nB)+nu+sum(G.nA));
    end;
   end;  % OK, now that derivative of f/resp wrt params available, get confidence regions
   GAMMA = [real(GAMMA);imag(GAMMA)];
   Q = GAMMA*P*GAMMA';  Gvar(k) = trace(Q);
   if (rem(k,kinc)==0)
    if det(Q) > 1e-20;  % Don't bother if error is negligible
     Qi = [0,-1;1,0]*Q*[0,1;-1,0]/det(Q);  % inv(Q) creates hassles if degenerate
     el = ellipse(Qi,level,pn); el = el(:);
     Ge(5:pn+4,k/kinc,r) = el + G.G(1,r,k)*ones(size(el));
     %  Now need to find max and min mag and phase when lying on these ellipses
     [v,l] = eig(Qi); v1 = v(:,1); v2 = v(:,2); l1 = l(1,1); l2 = l(2,2);
    else
     l1 = 1e20; l2 = l1;
    end;
    % Find lower and upper bounds on magnitude
    l=min(l1,l2); if (l>1e-10) l = sqrt(level/l); else l=1e10; end;
    minG = min(abs(G.G(1,r,:)));
    mupper(k/kinc) = abs(G.G(1,r,k))+l;
    mlower(k/kinc) = max(abs(G.G(1,r,k))-l,1e-2*minG);
    % That is bounds on mag done, now let's do bound on phase.
    if (l>abs(G.G(1,r,k)))
      per = pi;
    else
      per = asin(l/abs(G.G(1,r,k)));
    end;
    pupper(k/kinc) = aG(1,r,k)+per;
    plower(k/kinc) = aG(1,r,k)-per;
   end;
  end;

  % Pack error quantification data into output data structure
  Ge(1:4,1:floor(length(G.w)/kinc),r) = [mupper;mlower;pupper;plower];

  % Update index of block in theta pertaining to the G.B and G.A we are up to
  if ~strcmp(lower(G.type),'fir') Bind = Bind+G.nB(r)+1;  Aind = Aind+G.nA(r); end;
end; % END of loop over each input-output model










