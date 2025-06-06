%   SHOWBODE: This routine takes the results of an identification experiment
%   and presents the results as a Bode plots of the estimated system
%   together with estimated confidence regions (if they are available).
%
%   Usage is:
%
%   showbode(G,G1,G2...);
%
%   Where
%
%   G      = Data structure specifying estimated system as produced by
%            routines like est.m, barx.m, sid.m, onid.m, foe.m, farx.m etc.
%   G1...  = Stuctures, which must contain an element Gx.G whose
%            column(s) specify complex frequency responses that are to be
%            overlaid on the frequency response(s) in G. There must be as
%            many columns in each G1.G etc. as there are in G.G.  As
%            well, the structures G1... must contain elements Gx.w and
%            Gx.T which specify frequency axes and sampling periods.  The
%            typical case is that all these requirements are
%            automatically met by G1... being obtained by estimation
%            functions such as est.m, barx.m, foe.m etc that a that
%            automatically fill in these elements.
%   G.disp = Structure that can be used to customise the plotting style
%            according to the entries in:
%
% G.disp.colour    = Line colour
% G.disp.linestyle = Line style
% G.disp.legend    = Legend for plot
% G.disp.aux       = If set to 'magonly' or 'phonly' for the *first* plot then,
%                    respectively, only magnitude or phase plots are shown for
%                    first and subsequent plots
% G.disp.unit      = If set to 'rad' (default) then rad/s used as
%                    x-axis.  If set to 'hz' then Herz are used.
% G.disp.error     = If set nonzero (default) p5% CI error bounds are
%                    shown, but hidden if set to zero.
%
%
%   written by Brett Ninness,  School of EE & CS
%              Adrian Wills    University of Newcastle
%                              Australia.

% Copyright (C) Brett Ninness

function varargout = showbode(varargin)

%Set empty variable that will contain handles to figures
handle = [];

%Make sure we have an argument
lg=length(varargin);
if lg<1, error('Need at least one input argument'); end

%Call m2f for all models in case they don't already have a frequency response
for i=1:nargin,
	if ~isfield(varargin{i},'G')
  if ( isnumeric(varargin{i}) | isstruct(varargin{i}) )   
   varargin{i}=m2f(varargin{i});  % Compute freq response if given model structure or data
  else
   error('I cannot recognise what you are giving me - it is neither model or data?');
  end;
 end;
end
  
%Put G.G's into correct structure
for i=1:lg,
 x=size(varargin{i}.G);
 if length(x)<3,
  gnew=varargin{i}.G;
  varargin{i}.G=[];
  for k=1:size(gnew,2),
   varargin{i}.G(1,k,:)=gnew(:,k);
  end
 end
end

% Do some preliminary error checking
x=size(varargin{1}.G);
nu=x(2);
ny=1;
if length(x)<3 %  MIMO or MISO ?
 ny = 1;
else
 ny = x(1);
end

twopi=1; % Default if not otherwise set below
% What parameters have been given for primary Bode plot? Set rest to defaults.
if isfield(varargin{1},'disp'),
 if isfield(varargin{1}.disp,'aux'),
  if any(strcmp(varargin{1}.disp.aux,{'magonly','phonly'})),
   aux = varargin{1}.disp.aux;
  else
   error('I cannot recognise specification of G.disp.aux');
  end
 else
  aux='';
 end
 if isfield(varargin{1}.disp,'unit')
  if strcmpi(varargin{1}.disp.unit,'hz')
   twopi = 2*pi;
  elseif strcmpi(varargin{1}.disp.unit,'rad') twopi=1;
  else error('Unit specified for G.disp.unit is neither of rad or hz and hence not recognised');
  end;
 else
  twopi=1;
 end
else
 aux='';
end

%Set default colour and linestyle order
col=['b','r','k','g','m','c'];  % Set default set of colours
lin={'-','-.','--',':'};

% Define orientation of subplots
if (strcmp(aux,'magonly') | strcmp(aux,'phonly'))
 f1 = 111; f2 = 111;
else
 f1 = 211; f2 = 212;
end
%figure('Visible', 'off'); pindex = gcf;         % Find out the current figure open

%Start figures
for r=1:nu
 for k=1:ny
  % One figure per input-output pair
  %figure(pindex); handle = [handle pindex]; clf; pindex=pindex+1;
  if nargout == 1; handle = [handle figure('Visible', 'off', 'IntegerHandle','off')]; else figure; end;
  %Loop over number of models
  for i=1:lg,
   %Extract current model structure
   Gcur=varargin{i}; 
   
   %make sure Gcur.T is NOT equal to zero
   if isfield(Gcur,'T'),
    if Gcur.T<=0,
     Gcur.T=1;
    end
   else
    Gcur.T=1;
   end
   
   %Try to determine colour and linestyle
   if isfield(Gcur,'disp'),
    if isfield(Gcur.disp,'colour'),     c=Gcur.disp.colour;    else c = col(mod(i-1,length(col))+1); end
    if isfield(Gcur.disp,'linestyle'),  l=Gcur.disp.linestyle; else l = lin{mod(i-1,length(lin))+1}; end
   else
    c = col(mod(i-1,length(col))+1);
    l = lin{mod(i-1,length(lin))+1};
   end
   
   %Append legend information
   lgnd{i}='System response';
   if isfield(Gcur,'disp'),
    if isfield(Gcur.disp,'legend'),
     lgnd{i}=Gcur.disp.legend;
    end
   end
   
   if ~strcmp(aux,'phonly')
    %------------------------------------------------------------------
    % Magnitude Plot
    %------------------------------------------------------------------
    if(~isfield(Gcur, 'thumbnail'))
     subplot(f1);
    end
    if i>1, hold on; end
    semilogx(Gcur.w/twopi,20*log10(abs(squeeze(Gcur.G(k,r,:)))),[c,l],'linewidth',2);    
   end
   
   if ~strcmp(aux,'magonly')
    %------------------------------------------------------------------
    % Phase Plot
    %------------------------------------------------------------------
    subplot(f2);
    if i>1, hold on; end
    semilogx(Gcur.w/twopi,180/pi*unwrap(angle(squeeze(Gcur.G(k,r,:)))),[c,l],'linewidth',2);
   end
  end
  
  
  %------------------------------------------------------------------
  % Add in appropriate titles and axes labels and legend
  %------------------------------------------------------------------
  subplot(f1)
  if isfield(Gcur,'Ge')
   title(sprintf('Estimated system and confidence regions: input %i to output %i',r,k));
  else
   title(sprintf('Estimated system: input %i to output %i',r,k));
  end
  
  if ~strcmp(aux,'phonly'),
   subplot(f1); grid on;
   if twopi>1
    xlabel('Frequency (Hz)');
   else
    xlabel('Frequency (rad/s)');
   end;
   ylabel('Magnitude (dB)');
   legend(lgnd,'Location','southwest');      
  end
  if ~strcmp(aux,'magonly'),
   subplot(f2); grid on;
   if twopi>1
    xlabel('Frequency (Hz)');
   else
    xlabel('Frequency (rad/s)');
   end;
   ylabel('Angle (deg)');
   legend(lgnd,'Location','southwest');      
  end
  
  
  %------------------------------------------------------------------
  % Finally, add error bound information if it is available
  %------------------------------------------------------------------
  if isfield(Gcur,'type'),
   if any(strcmpi(Gcur.type,{'fir','arx','armax','bj','oe'})),
    err=0;
    if isfield(Gcur,'disp'),
     if isfield(Gcur.disp,'error'),
      err = Gcur.disp.error;
     end
    end
    if err,
     try
      if ~isfield(Gcur,'Ge')   % Don't over-write! Maybe MCMC already filled this bit in
       [Gcur.Ge, Gcur.var] = gerror(Gcur);
      end;
      
      % Extract upper and lower magnitude and phase bounds
      mupper = Gcur.Ge(1,:,r); mlower = Gcur.Ge(2,:,r);
      pupper = Gcur.Ge(3,:,r); plower = Gcur.Ge(4,:,r);
      
      %  The bounds are provided at subsampled points of the estimated
      %  frequency respose G - so get a frequency axis that reflects this
      %  subsampling
      kinc = floor(length(Gcur.w)/length(mupper));
      kk = kinc:kinc:kinc*floor(length(Gcur.w)/kinc);
      ww = Gcur.w(kk);
      % Plot of magnitude error bounds
      if norm(mupper-mlower)>1e-6,
       if ~strcmp(aux,'phonly')
        subplot(f1);
        hold on;
        semilogx(ww,20*log10([mupper(:),mlower(:)]),'-.b');
        hf=fill([ww(:);flipud(ww(:))],20*log10([mupper(:);flipud(mlower(:))]),[183,209,209]/255);
        alpha(0.01);
        hold off;
       end;
      end
      % Plot phase error bounds
      if norm(pupper-plower)>1e-6,
       if ~strcmp(aux,'magonly')
        subplot(f2);
        hold on;
        semilogx(ww,180/pi*[pupper(:),plower(:)],'-.b');
        hf=fill([ww(:);flipud(ww(:))],180/pi*[pupper(:);flipud(plower(:))],[183,209,209]/255);
        alpha(0.01);
        hold off;
       end;
      end
     catch
      disp('Could not display error bounds');
     end
    end
   end;
  end
 end; % End of loop over outputs
end; % End of loop over inputs

%Make sure we map the handle to the output if needed
if nargout==1,
 varargout{1} = handle;
end

hold off;  % Just to make sure no plots being held.
