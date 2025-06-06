%   NONPAR This routine takes a record of input-output data [y,u] and
%   returns a non-parametric estimate of the frequency response G(e^{j*w})
%   of the linear system the might have produced this data according to
%
%   y_t = G(q)u_t
%
%   Usage is:
%
%   G = nonpar(Z,M,OPT);
%
%   where
%
%   Z         = Input-Output data in one of two forms.  The standard form
%               is for it to be a record with elements Z.y and Z.u, each
%               of which are matrices with number of rows equal to the
%               number of data samples, and number of columns equal (respectively)
%               to the number of outputs and the number of inputs.  On
%               the other hand, Z can be a matrix of the form Z = [y,u]
%               where it is assumed that y is a column vector of output
%               measurements and u is a matrix whose columns are the
%               input measurements; in this latter MISO models are
%               being considered.
%   M         = Data structure which defines the model structure which
%               is to be estimated from the data as follows:
%    M.w      = vector of frequencies at which to calculate frequency
%               response of estimated model.  Specify in real frequency,
%               not normalised.  Default is 3 decades up to folding freq.
%    M.delay  = Number of samples of delay to include. In the
%               case of a MIMO system, this should be a vector of delays,
%               one for each input being considered.
%   OPT       = Data structure which defines options for the estimation
%               algorithm as follows:
%    OPT.n    = number of starting data points to discard to get
%               rid of initial condition effects.  Default is none.
%    OPT.alg  = The algorithm type used.  It may be set as:
%
% 	      'bltuk' - "Blackman-Tukey" method in which the estimate is
% 	                 as the ratio \Phi_{yu}(w)/\Phi_u(w) of
% 	                 cross-spectrum estimate to input spectrum
% 	                 estimate, and these are found as DFT's of
% 	                 estimated cross-covariance and covariance.
%                        This is the default algorithm.
% 	      'etfe'   - "Empirical Transfer Function Estimate" in which
% 	                 result is found as ratio of input output DFT to input DFT.
%    OPT.window    - Window function used by either of above
%                        methods. It may be set to any of `boxcar',
%                        `bartlett',  `hamming' or `hanning'.
%                        Default is `hanning'.
%    OPT.N         - Window length used in Blackman Tukey method.
%                    The default is 20% of the data length.
%
%   G         = Data structure which specifies the estimated model as
%               follows:
%    G.G      = Estimated Non-parametric frequency response.
%
%    Written by Brett Ninness, School of EE & CS
%                              University of Newcastle
%                          Australia.

% Copyright (C) Brett Ninness.

function G = nonpar(Z,M,OPT);

%Call startZ
Z = startZ(Z);

% Extract input and output from data matrix
[y,u,ny,nu,Ny] = Z2data(Z);

% Unspecified parts of M -> defaults
M.type = 'nonpar'; 
M      = startM(Z,M);
M.w    = M.w*M.T;
if isfield(Z,'disp') && ~isfield(M,'disp'),
 M.disp = Z.disp; 
end
G      = M;

%Sample time must come from the data
M.T = Z.T;

%Switch between frequency data and time domain data
switch Z.type,
 case 'frequency'
  G.G = Z.y;
  G.w = Z.w;
  G.disp.legend = 'Raw Data';
  
 case 'time'
  % Include delays specified in model structure on inputs
  for r=1:nu u(:,r) = [zeros(M.delay(r),1);u(1:Ny-M.delay(r),r)]; end;
  
  % Force data lengths to be even
  Ny2 = floor(Ny/2); Ny = 2*Ny2; y = y(1:Ny); if nu>0 u = u(1:Ny); end;
  
  % Unspecified parts of OPT -> defaults
  if ~exist('OPT') OPT = startOPT([]); else OPT = startOPT(OPT);        end;
  if (OPT.n>=Ny) error('Cannot OPT.n larger than height of Z!');        end;
  if ~isfield(OPT,'alg')         OPT.alg    = 'bltuk';                  end;
  if strcmpi(OPT.alg,'gn')       OPT.alg    = 'bltuk';                  end;
  if ~isfield(OPT,'window')      OPT.window = 'hanning';                end;
  if ~isfield(OPT,'N')           OPT.N      =  max(1,floor(0.2*Ny));    end;
  
  
  % Force data windowing length to be odd
  Nwin2 = floor(OPT.N/2); OPT.N = 2*Nwin2+1;
  
  % Calculate preliminaries for windowing sequence to be used on data
  k = -Nwin2:1:Nwin2; wun = ones(size(k)); side = zeros(1,Ny2-Nwin2);
  
  
  % Check for type of algorithm specified
  switch lower(OPT.alg)
   
   case {'etfe'},  % Empirical Transfer Function Estimate Selected
    
    lbl = 'ETFE';
    
    switch lower(OPT.window)
     case {'hamming'},  win = [side,0.54*wun+0.46*cos(pi*k./OPT.N) side];
     case {'hanning'},  win = [side,0.5*wun+0.5*cos(pi*k./OPT.N) side];
     case {'bartlett'}, win = [side,(OPT.N-abs(k))/OPT.N,side];
     case {'boxcar'},   win = ones(1,Ny);
     otherwise
      error('What sort of window is that?  Should be one of "hanning", "hamming", "bartlett" or "boxcar"');
    end;
    win = win(1:Ny); win = win(:);
    
    % Window the data before calculating spectrums - same as smoothing ETFE afterwards
    yw = y.*win; if nu>0 uw = u.*win; end;
    % Form non-parametric estimate as ratio of input and output spectra
    Y = dft(yw,M.w); if nu>0 U = dft(uw,M.w); G.G = Y./U; else G.G = Y; end;
    G.G(isnan(G.G)) = 0;  % Check for divide by zero.
    
    
   case {'bltuk'}  % Blackman-Tukey method selected
    
    lbl = 'Blackman-Tukey';
    
    % Calculate windowing sequence to be used on data
    switch lower(OPT.window)
     case {'hanning'},  win = 0.5*wun +0.5*cos(pi*k./Nwin2);
     case {'hamming'},  win = 0.54*wun+0.46*cos(pi*k./Nwin2);
     case {'bartlett'}, win = (Nwin2*wun-abs(k))/Nwin2;
     case {'boxcar'},   win = ones(1,2*Nwin2+1);
     otherwise error('What sort of window is that?  Should be one of "hanning", "hamming", "bartlett" or "boxcar"');
    end;    
    
    % Calculate input covariance sequence and input-output cross covariance
    Y = fft(y); Y = Y(:);
    if nu>0
     U = fft(u); U = U(:);
     Ru = real(ifft(abs(U).^2))/Ny; Ryu = real(ifft(Y.*conj(U)))/Ny;
     
     % Because FFT uses 0->2pi rather than -pi->pi, then neg lag correls are at end.  Fix.
     RRu = [Ru(Ny-Nwin2+1:Ny);Ru(1:Nwin2+1)];
     RRyu = [Ryu(Ny-Nwin2+1:Ny);Ryu(1:Nwin2+1)];
     G.G = dft(RRyu.*win(:),M.w)./dft(RRu.*win(:),M.w);
    else
     Ry = real(ifft(abs(Y).^2))/Ny; RRy = [Ry(Ny-Nwin2+1:Ny);Ry(1:Nwin2+1)];
     G.G = dft(RRy.*win(:),M.w); tmp = exp(-j*M.w*Nwin2); G.G = G.G(:).*tmp(:);
    end;
    
   otherwise
    error('What sort of estimation method is that? OPT.alg should be one of "bltuk" or "etfe"')
    
  end;
  
  % Include effect of any delays specified in model structure
  pdel = exp((-j*M.w(:)*M.T)*M.delay');
  G.G = G.G(:).*pdel(:,1);
  G.w = M.w(:)/Z.T; G.delay = M.delay; G.T = M.T;
  G.disp.legend = ['Nonparametric ',lbl,' Estimate'];
end


% Set some final arguments for compatibility
G.var = 0;
G.alg = lbl;
G.OPT = OPT;
G.Ny  = Ny;


















