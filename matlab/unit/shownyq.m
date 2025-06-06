%   SHOWNYQ: This routine takes the results of an identification experiment and
%   presents the results as a Nyquist plot(s) of the estimated system(s)
%   together with estimated confidence regions (if they are available).
%
%   Usage is:
%
%   shownyq(G,G1,G2...);
%
%   Where
%
%   G     =  Data structure specifying estimated system as produced by
%            routines like est.m, barx.m, onid.m.
%   G1... =  Stuctures, which must contain an element Gx.G whose
%            column(s) specify complex frequency responses that are to be
%            overlaid on the frequency response(s) in G. There must be as
%            many columns in each G1.G etc. as there are in G.G.  As
%            well, the structures G1... must contain elements Gx.w and
%            Gx.T which specify frequency axes and sampling periods.  The
%            typical case is that all these requirements are
%            automatically met by G1..G8 being obtained by estimation
%            functions such as est.m, barx.m, foe.m etc that a that
%            automatically fill in these elements.
%   G.disp = Structure that can be used to customise the plotting style
%            according to the entries in"
%
% G.disp.colour    = Line colour
% G.disp.linestyle = Line style
% G.disp.legend    = Legend for plot
% G.disp.axes      = Axes handle for plot
% G.disp.error     = If set nonzero (default) p5% CI error bounds are
%                    shown, but hidden if set to zero.
%
%
%   written by Brett Ninness  School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness

function handle = shownyq(varargin)

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
  gnew=varargin{i}.G(:);
  varargin{i}.G=[];
  varargin{i}.G(1,1,:)=gnew;
 end
end

%keyboard

% Do some preliminary error checking
if ~isstruct(varargin{1}) error('Input model must be a structure and not just a matrix'); end;
x=size(varargin{1}.G); nu=x(2); ny=1; if length(x)<3 ny=1; else ny = x(1); end; %  MIMO or MISO ?
if ~isfield(varargin{1},'G') error('Input must be in model structure form'); end;

%Set default colour and linestyle order
col=['b','r','k','g','m','c'];  % Set default set of colours
lin={'-','-.','--',':'};

%figure; pidx=gcf;
figure('Visible', 'off', 'IntegerHandle','off'); pidx=gcf;
for r=1:nu  % Loop through all inputs
 %figure(r+pidx-1); handle = [handle r+pidx-1]; clf;           % One figure per input
 if nargout == 1; handle = [handle figure('Visible', 'off', 'IntegerHandle','off')]; else figure; end;
 for k=1:ny
  subplot(ny,1,k);
  
  %Loop over number of models
  for i=1:lg,
   %Extract current model structure
   Gcur=varargin{i};
   
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
   
   %--------------------------------------------------------------
   % Nyquist Plot
   %--------------------------------------------------------------
   if i>1, hold on; end
   plot(real(squeeze(Gcur.G(k,r,:))),imag(squeeze(Gcur.G(k,r,:))),[c,l],'linewidth',2);
   
   
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
       [m,dummy] = size(Gcur.Ge(:,:,r));
       ell = Gcur.Ge(5:m,:,r);
       hold on
       plot(real(ell),imag(ell),'-.b');
       hf=fill(real(ell),imag(ell),[183,209,209]/255);
       alpha(0.5);
       hold off
      catch
       disp('Could not display error ellipsoids')
      end
     end
    end
   end
  end
  
  
  %------------------------------------------------------------------
  % Add in appropriate titles and axes labels and legend
  %------------------------------------------------------------------
  if isfield(varargin{1},'Ge')
   title(sprintf('Estimated system and confidence regions: input %i to output %i',r,k));
  else
   title(sprintf('Estimated system: input %i to output %i',r,k));
  end
  legend(lgnd,'Location','northeast');   
  grid on;
  xlabel('Real');
  ylabel('Imaginary');
  
 end; % Loop over outputs
end; % Loop over inputs

hold off;  % Just to make sure no plots being held.
