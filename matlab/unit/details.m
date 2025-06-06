% DETAILS: This function summarises details about a model structure.
%
% Usage is
%
% details(G)
%
% Where
%
% G = data structure storing an estimated model.
%
% Example: G = est(Z,M,OPT); details(G);
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness.

function handle = details(varargin)

%Set empty variable that will contain handles to figures
handle = [];
popup  = [];

lg=length(varargin);
if lg<1, error('Need at least one input argument'); end
G=varargin{1};
if nargin>1,
 popup = varargin{2};
end

%Detect if gui is running
gui = 0;
guih = [];
if isfield(G,'OPT')   % Hack by Brett - should really fix other routines so G.OPT is always put in
 if isfield(G.OPT,'gui'),
  if ~isempty(G.OPT.gui)
   gui  = 1;           %GUI is running
   guih = G.OPT.gui;   %GUI handle
  end
 end
else
 gui=0;
end

urline = '----------------------------------------------------------------------';
urdisp(' ',gui,guih,popup);

switch G.type
 
 %----------------------------------------------------------------------
 %                                NONPAR
 %----------------------------------------------------------------------
 case 'nonpar'
  %Try a quick return if NONPAR
  if strcmpi(G.type,'nonpar'),
   urdisp(urline,gui,guih,popup);
   urdisp(['Details for Estimated Model Structure'],gui,guih,popup);
   urdisp(urline,gui,guih,popup);
   urdisp(['Model Structure Used           = Non Parametric'],gui,guih,popup)
   if strcmpi(G.alg,'ETFE')
    urdisp(['Estimation algorithm           = ETFE'],gui,guih,popup)
   else
    urdisp(['Estimation algorithm           = Blackman-Tukey'],gui,guih,popup)
   end
   urdisp(sprintf('Data length (samples)          = %i',G.Ny),gui,guih,popup)
   urdisp(sprintf('Window used                    = %s',G.OPT.window),gui,guih,popup)
   urdisp(sprintf('Window length (samples)        = %i  (%3.1f%%)',G.OPT.N,100*G.OPT.N/G.Ny),gui,guih,popup)
   urdisp(urline,gui,guih,popup)
   return;
  end
  
  
  %----------------------------------------------------------------------
  %                      TRANSFER FUNCTION MODELS
  %----------------------------------------------------------------------
 case {'ar','arma','fir','arx','armax','oe','bj'}
  
  if ~isfield(G,'var')
   error('This function is to be used with estimated models, not just structure specifications!')
  else  % Otherwise, get standard deviations on estimated parameters
   % Need to discriminate several times on type - streamline the syntax
   if strcmpi(G.type,'armax') armax=1; else armax=0; end;
   if strcmpi(G.type,'fir')   fir=1;   else fir=0;   end;
   
   % Need to handle fir/arx as specials since latter hand back theta/SD with no nonlin parts
   % Are there any non-linear parts on inputs?
   inlin=1;
   for k=1:length(G.in)
    if ~strcmp(lower(G.in(k).type),'linear') inlin=0; end
   end
   if any(strcmpi(G.type,{'arx','farx','fir'}))
    SD = G.SD; % SD's of linear parameter parts
    if ~inlin
     th  = m2theta(G);
     g   = G;
     g.P = g.P(end-length(th)+1:end,end-length(th)+1:end);
     [dummy,sd] = theta2m(th,g);
     SD.in = sd.in;
    end
    if ~strcmpi(G.out.type,'linear')
     th  = m2theta(G);
     g   = G;
     g.P = g.P(end-length(th)+1:end,end-length(th)+1:end);
     [dummy,sd] = theta2m(th,g);
     SD.out = sd.out;
    end
   else % Otherwise, all parameters found iteratively, so theta2m extracts all SD's
    [dummy,SD] = theta2m(G.theta,G);
   end
  end
  
  if ~isfield(G,'type'),
   G.type = 'not specified';
  end
  
  [nu,dummy,ny] = size(G.A);  % How many dynamics models to print out?
  
  %  Print out the details of the linear part of the model structure:
  
  urdisp(urline,gui,guih,popup);
  urdisp(['Details for Estimated Model Structure'],gui,guih,popup);
  urdisp(urline,gui,guih,popup);
  urdisp(['Operator used in model         = ',G.op],gui,guih,popup)
  if G.T>0,
   urdisp(['Sampling Period                = ',sprintf('%f ',G.T),'seconds'],gui,guih,popup)
  else
   urdisp(['Sampling Period                = irregular'],gui,guih,popup)
  end
  if G.T>0,
   urdisp(['Estimated Innovations Variance = ',sprintf('%e ',trace(G.var))],gui,guih,popup)
  else
   urdisp(['Estimated Innovations Variance = ',sprintf('%e ',trace(G.ss.R))],gui,guih,popup)
  end
  switch lower(G.type),
   case {'fir','nfir'}
    urdisp(['Model Structure Used           = FIR'],gui,guih,popup)
   case {'arx','narx'}
    urdisp(['Model Structure Used           = ARX'],gui,guih,popup);
   case {'bj','nbj'}
    urdisp(['Model Structure Used           = Box-Jenkins'],gui,guih,popup)
   case {'arma','narma'}
    urdisp(['Model Structure Used           = ARMA'],gui,guih,popup)
   case {'armax','narmax'}
    urdisp(['Model Structure Used           = ARMAX'],gui,guih,popup)
   case {'oe','noe'}
    urdisp(['Model Structure Used           = Output Error'],gui,guih,popup)
   otherwise
    urdisp(['Model Structure Used           = Not Specified'],gui,guih,popup)
  end
  
  if strcmpi(G.alg,'block')
   urdisp(['Estimation algorithm           = Block solution'],gui,guih,popup)
  elseif strcmpi(G.alg,'gn')
   urdisp(['Estimation algorithm           = Gauss-Newton search'],gui,guih,popup);
  elseif strcmpi(G.alg,'em')
   urdisp(['Estimation algorithm           = Expect-Maxim (EM) search'],gui,guih,popup)
  elseif strcmpi(G.alg,'sid')
   urdisp(['Estimation algorithm           = N4SID Subspace Regression'],gui,guih,popup)
  elseif strcmpi(G.alg,'cca')
   urdisp(['Estimation algorithm           = Canonical Correlation Analysis'],gui,guih,popup)
  elseif strcmpi(G.alg,'n4sid')
   urdisp(['Estimation algorithm           = N4SID Subspace Regression'],gui,guih,popup)
  else
   urdisp(['Estimation algorithm           = Not Specified'],gui,guih,popup)
  end;
  
  % Print out details for non-linear part of model structure
  if G.nu>0  % Only if there is an input
   for k=1:nu
    urdisp([sprintf('Input #%d block type            = ',k),G.in(k).type],gui,guih,popup)
   end;
   urdisp(['Output block type              = ',G.out.type],gui,guih,popup)
   urdisp(urline,gui,guih,popup)
  end;
  
  %Print the model equations
  if isfield(G,'modelEquations'),
   urdisp(sprintf('\nModel Equations are:'),gui,guih,popup)
   urdisp(G.modelEquations,gui,guih,popup)
  end
  
  if ny*nu < 20  % Don't print out a ridiculous number of transfer functions
   for m=1:ny  % Loop over all outputs
    for k=1:nu  % Loop over all inputs
     urdisp(' ',gui,guih,popup)
%      if fir
%       urdisp(urline,gui,guih,popup);
%       urdisp(['Estimated (generalised) FIR parameters & standard deviations:'],gui,guih,popup);
%       urdisp(urline,gui,guih,popup);
%       th='TH= ';sd='SD= ';
%       for r=1:length(G.TH(:,k));
%        %  theta = [G.TH(:,k)';SD.th(:,k)']
%        if (G.TH(r,k)<0) sp1=''; else sp1=' '; end;
%        th = [th,sp1,sprintf('%0.4f ',G.TH(r,k))];
%        sd = [sd,sprintf(' %0.4f ',SD.th(r,k))];
%       end;
%       urdisp(th,gui,guih,popup); urdisp(sd,gui,guih,popup); urdisp(' ',gui,guih,popup);
%      end;
     
     urdisp(urline,gui,guih,popup);
     urdisp(sprintf('Input #%d to Output #%d Estimated T/F model + standard devs: ',k,m),gui,guih,popup);
     urdisp(urline,gui,guih,popup);
     
     if strcmpi(G.type,'arma'),
      B1=[    ]; B2='C = ';B3='SD= ';
      SD.C = [0 SD.C];
      for r=1:G.nC(k)+1
       if r>1
        B1 = [B1,[' ',G.op,'^-',int2str(r-1),'   ']];
       else
        B1 = [B1,['     1   ','   ']];
       end;
       if (G.C(r)<0) sp1=''; else sp1=' '; end;
       B2 = [B2,sp1,sprintf('%0.4f ',G.C(r))];
       B3 = [B3,sprintf(' %0.4f ',SD.C(r))];
      end;
      urdisp(B1,gui,guih,popup);
      urdisp(B2,gui,guih,popup);
      urdisp(B3,gui,guih,popup);
     else
      B1=[    ]; B2='B = ';B3='SD= ';
      if ~isempty(G.B) % Could be here in AR case => no B to print
       for r=1:G.nB(k)+1
        if r>1
         if strcmpi(G.op,'s')
          B1 = [B1,[' ',G.op,'^',int2str(G.nB(k)-r+1),'   ']];
         else 
          B1 = [B1,[' ',G.op,'^-',int2str(r-1),'   ']];
         end;
        else
         if strcmpi(G.op,'s')
          B1 = [B1,['     ',G.op,'^',int2str(G.nB(k)),'   ',' ']];
         else
          B1 = [B1,['     1   ','   ']];
         end;
        end;
        if (G.B(k,r)<0) sp1=''; else sp1=' '; end;
        if strcmpi(G.op,'s')
         dB = G.nA(k)-G.nB(k);
         B2 = [B2,sp1,sprintf('%2.3f ',G.B(k,dB+r,m))];
        else
         B2 = [B2,sp1,sprintf('%0.4f ',G.B(k,r,m))];         
        end;
        if [~fir] B3 = [B3,sprintf(' %0.4f ',SD.B(k,r))]; end;
       end;
       urdisp(B1,gui,guih,popup); urdisp(B2,gui,guih,popup); if [~fir] urdisp(B3,gui,guih,popup); end
      end
     end
     A1=[    ]; A2='A = ';A3='SD=  0      ';
     for r=1:G.nA(k)+1
      if r>1
       if strcmpi(G.op,'s')       
        A1 = [A1,[' ',G.op,'^',int2str(G.nA(k)-r+1),'   ']];
       else
        A1 = [A1,[' ',G.op,'^-',int2str(r-1),'   ']];        
       end;
      else
       if strcmpi(G.op,'s')           
        A1 = [A1,['     ',G.op,'^',int2str(G.nA(k)),'   ',' ']];        
       else 
        A1 = [A1,['     1   ','   ']];
       end;        
      end;
      if (G.A(k,r)<0) sp1=''; else sp1=' '; end;
      if strcmpi(G.op,'s')
       A2 = [A2,sp1,sprintf('%2.3f ',G.A(k,r,m))];       
      else
       A2 = [A2,sp1,sprintf('%0.4f ',G.A(k,r,m))];
      end;       
      if [r>1,~fir] A3 = [A3,sprintf(' %0.4f ',SD.A(k,r-1))]; end;
     end;
     urdisp(' ',gui,guih,popup);
     urdisp(A1,gui,guih,popup); urdisp(A2,gui,guih,popup); if [~fir] urdisp(A3,gui,guih,popup); end; urdisp(' ',gui,guih,popup);
     urdisp(['delay = ',int2str(G.delay(k)), ' samples'],gui,guih,popup); urdisp(' ',gui,guih,popup);
     
     % Print out where the poles are
     p = cplxpair(roots(G.A(k,1:G.nA(k)+1)));
     idx = 1; pls = [];
     while idx<=G.nA(k)
      mg = abs(p(idx)); ph = abs(angle(p(idx)));  % Complex conj pair?
      if abs(imag(p(idx)))>0
       idx = idx+2;
       if idx<=G.nA(k) cm=','; else cm='.'; end;
       pls = [pls,' ',sprintf('%0.4f*exp(+-j%0.4f)',mg,ph),cm];
      else
       idx = idx+1;
       if idx<=G.nA(k) cm=','; else cm='.'; end;
       pls = [pls,' ',sprintf('%0.4f',p(idx-1)),cm];
      end;
     end;
     urdisp(['Poles at',pls],gui,guih,popup)
     %    urdisp(' ');
    end;  % End of loop over inputs
   end; % End of loop over outputs
  else % If too many input output pairs
   urdisp(urline,gui,guih,popup);
   urdisp(sprintf('%d Inputs & %d Outputs => %d T/Fs: Too many to urdisplay! ',nu,ny,nu*ny),gui,guih,popup);
   urdisp(urline,gui,guih,popup);
  end;
  
  if any(strcmpi(G.type,{'bj','armax'}))
   urdisp(' ',gui,guih,popup);
   for m=1:ny
    for k=1:ny
     urdisp(urline,gui,guih,popup);
     urdisp(sprintf('Innovations Input #%d to Output #%d Estimated Spec Fact + stand dev: ',k,m),gui,guih,popup);
     urdisp(urline,gui,guih,popup);
     
     C1=[    ]; C2='C = ';C3='SD=  0      ';   D2='D = ';D3=C3;
     for r=1:length(G.C(m,:,k))
      if r>1
       C1 = [C1,[' ',G.op,'^-',int2str(r-1),'   ']];
      else
       C1 = [C1,['     1   ','   ']];
      end;
      if (G.C(1,r)<0) sp1=''; else sp1=' '; end;
      C2 = [C2,sp1,sprintf('%0.4f ',G.C(1,r))];
      if [r>1] C3 = [C3,sprintf(' %0.4f ',SD.C(1,r-1))]; end;
      if ~armax
       if (G.D(1,r)<0) sp2=''; else sp2=' '; end;
       D2 = [D2,sp2,sprintf('%0.4f ',G.D(1,r))];
       if [r>1] D3 = [D3,sprintf(' %0.4f ',SD.D(1,r-1))]; end;
      end;
     end;
     
     urdisp(' ',gui,guih,popup); urdisp(C1,gui,guih,popup); urdisp(C2,gui,guih,popup);
     urdisp(C3,gui,guih,popup);
     urdisp(' ',gui,guih,popup);
     if ~armax
      urdisp(' ',gui,guih,popup);
      urdisp(C1,gui,guih,popup);
      urdisp(D2,gui,guih,popup);
      urdisp(D3,gui,guih,popup);
      urdisp(' ',gui,guih,popup);
     end
     
     % Print out where the poles/zeros of the noise model are - first for C:
     p = cplxpair(roots(G.C(m,:,1))); idx = 1; pls = []; np = length(p);
     while idx<=np
      mg = abs(p(idx)); ph = abs(angle(p(idx)));  % Complex conj pair?
      if abs(imag(p(idx)))>0
       idx = idx+2; if idx<=np cm=','; else cm='.'; end;
       pls = [pls,' ',sprintf('%0.4f*exp(+-j%0.4f)',mg,ph),cm];
      else
       idx = idx+1; if idx<=np cm=','; else cm='.'; end;
       pls = [pls,' ',sprintf('%0.4f',p(idx-1)),cm];
      end;
     end;
     urdisp(['Zeros of C at:',pls],gui,guih,popup)
     urdisp(' ',gui,guih,popup);
     
     % Now urdisplay zeros of D
     p = cplxpair(roots(G.D(m,:,1))); idx = 1; pls = []; np = length(p);
     while idx<=np
      mg = abs(p(idx)); ph = abs(angle(p(idx)));  % Complex conj pair?
      if abs(imag(p(idx)))>0
       idx = idx+2; if idx<=np cm=','; else cm='.'; end;
       pls = [pls,' ',sprintf('%0.4f*exp(+-j%0.4f)',mg,ph),cm];
      else
       idx = idx+1; if idx<=np cm=','; else cm='.'; end;
       pls = [pls,' ',sprintf('%0.4f',p(idx-1)),cm];
      end;
     end;
     urdisp(['Zeros of D at:',pls],gui,guih,popup)
     urdisp(' ',gui,guih,popup);
    end; % End of loop over input innovation
   end;  % End of loop over output number
  end;  % End of check on whether there is a noise model worth outputting
  
  
  
  %----------------------------------------------------------------------
  %                            STATE SPACE MODELS
  %----------------------------------------------------------------------
 case {'ss','bilin','bilinear'}
  urdisp(urline,gui,guih,popup);
  urdisp(['Details for Estimated Model Structure'],gui,guih,popup);
  urdisp(urline,gui,guih,popup);
  urdisp(['Operator used in model         = ',G.op],gui,guih,popup)
  if G.T>0,
   urdisp(['Sampling Period                = ',sprintf('%f ',G.T),'seconds'],gui,guih,popup)
  else
   urdisp(['Sampling Period                = irregular'],gui,guih,popup)
  end
  if G.T>0,
   urdisp(['Estimated Innovations Variance = ',sprintf('%e ',trace(G.var))],gui,guih,popup)
  else
   urdisp(['Estimated Innovations Variance = ',sprintf('%e ',trace(G.ss.R))],gui,guih,popup)
  end
  switch lower(G.type),
   case {'ss','nss'}
    urdisp(['Model Structure Used           = Linear State Space'],gui,guih,popup)
   case {'bilin','bilinear'}
    urdisp(['Model Structure Used           = Bilinear State Space'],gui,guih,popup)
   otherwise
    urdisp(['Model Structure Used           = Not Specified'],gui,guih,popup)
  end
  
  if strcmpi(G.alg,'gn')
   urdisp(['Estimation algorithm           = Gauss-Newton search'],gui,guih,popup);
  elseif strcmpi(G.alg,'em')
   urdisp(['Estimation algorithm           = Expect-Maxim (EM) search'],gui,guih,popup)
  elseif strcmpi(G.alg,'sid')
   urdisp(['Estimation algorithm           = N4SID Subspace Regression'],gui,guih,popup)
  elseif strcmpi(G.alg,'cca')
   urdisp(['Estimation algorithm           = Canonical Correlation Analysis'],gui,guih,popup)
  elseif strcmpi(G.alg,'n4sid')
   urdisp(['Estimation algorithm           = N4SID Subspace Regression'],gui,guih,popup)
  else
   urdisp(['Estimation algorithm           = Not Specified'],gui,guih,popup)
  end;
  
  %Print the model equations
  if isfield(G,'modelEquations'),
   urdisp(sprintf('\nModel Equations are:'),gui,guih,popup)
   urdisp(G.modelEquations,gui,guih,popup)
  end
  
  urdisp(urline,gui,guih,popup);
  urdisp('    State Space Matrices ',gui,guih,popup);
  urdisp(urline,gui,guih,popup);
  
  frmt = '%11.2e';
  if numel(G.ss.A)>0,
   str = num2str(G.ss.A,frmt);
   if any(strfind(str(:,1)','-')),
    lbl1 = 'A  = ';
    lbl2 = '     ';
   else
    lbl1 = 'A  =  ';
    lbl2 = '      ';
   end
   urdisp([lbl1 str(1,:)],gui,guih,popup);
   for i=2:size(str,1)
    urdisp([lbl2 str(i,:)],gui,guih,popup);
   end
   urdisp(' ',gui,guih,popup);
  end
  if numel(G.ss.B)>0,
   str = num2str(G.ss.B,frmt);
   if any(strfind(str(:,1)','-')),
    lbl1 = 'B  = ';
    lbl2 = '     ';
   else
    lbl1 = 'B  =  ';
    lbl2 = '      ';
   end
   urdisp([lbl1 str(1,:)],gui,guih,popup);
   for i=2:size(str,1)
    urdisp([lbl2 str(i,:)],gui,guih,popup);
   end
   urdisp(' ',gui,guih,popup);
  end
  if numel(G.ss.C)>0,
   str = num2str(G.ss.C,frmt);
   if any(strfind(str(:,1)','-')),
    lbl1 = 'C  = ';
    lbl2 = '     ';
   else
    lbl1 = 'C  =  ';
    lbl2 = '      ';
   end
   urdisp([lbl1 str(1,:)],gui,guih,popup);
   for i=2:size(str,1)
    urdisp([lbl2 str(i,:)],gui,guih,popup);
   end
   urdisp(' ',gui,guih,popup);
  end
  if numel(G.ss.D)>0,
   str = num2str(G.ss.D,frmt);
   if any(strfind(str(:,1)','-')),
    lbl1 = 'D  = ';
    lbl2 = '     ';
   else
    lbl1 = 'D  =  ';
    lbl2 = '      ';
   end
   urdisp([lbl1 str(1,:)],gui,guih,popup);
   for i=2:size(str,1)
    urdisp([lbl2 str(i,:)],gui,guih,popup);
   end
   urdisp(' ',gui,guih,popup);
  end
  if numel(G.ss.K)>0,
   str = num2str(G.ss.K,frmt);
   if any(strfind(str(:,1)','-')),
    lbl1 = 'K  = ';
    lbl2 = '     ';
   else
    lbl1 = 'K  =  ';
    lbl2 = '      ';
   end
   urdisp([lbl1 str(1,:)],gui,guih,popup);
   for i=2:size(str,1)
    urdisp([lbl2 str(i,:)],gui,guih,popup);
   end
   urdisp(' ',gui,guih,popup);
  end
  if numel(G.ss.F)>0,
   str = num2str(G.ss.F,frmt);
   if any(strfind(str(:,1)','-')),
    lbl1 = 'F  = ';
    lbl2 = '     ';
   else
    lbl1 = 'F  =  ';
    lbl2 = '      ';
   end
   urdisp([lbl1 str(1,:)],gui,guih,popup);
   for i=2:size(str,1)
    urdisp([lbl2 str(i,:)],gui,guih,popup);
   end
   urdisp(' ',gui,guih,popup);
  end
  if numel(G.ss.G)>0,
   str = num2str(G.ss.G,frmt);
   if any(strfind(str(:,1)','-')),
    lbl1 = 'G  = ';
    lbl2 = '     ';
   else
    lbl1 = 'G  =  ';
    lbl2 = '      ';
   end
   urdisp([lbl1 str(1,:)],gui,guih,popup);
   for i=2:size(str,1)
    urdisp([lbl2 str(i,:)],gui,guih,popup);
   end
   urdisp(' ',gui,guih,popup);
  end
  if numel(G.ss.X1)>0,
   str = num2str(G.ss.X1,frmt);
   if any(strfind(str(:,1)','-')),
    lbl1 = 'X1 = ';
    lbl2 = '     ';
   else
    lbl1 = 'X1 =  ';
    lbl2 = '      ';
   end
   urdisp([lbl1 str(1,:)],gui,guih,popup);
   for i=2:size(str,1)
    urdisp([lbl2 str(i,:)],gui,guih,popup);
   end
   urdisp(' ',gui,guih,popup);
  end
  if numel(G.ss.A)>0,
   pls = eig(G.ss.A)';
   if isreal(pls),
    urdisp(['Poles of A  = ' num2str(pls,frmt)],gui,guih,popup);
   else
    %urdisp(['Poles of A             = ' num2str(pls,frmt)],gui,guih,popup);
    urdisp(['Poles of A (magnitude) = ' num2str(abs(pls),frmt)],gui,guih,popup);
   end
  end
end

%Extract number of inputs and outputs
nu = G.nu;
ny = G.ny;


%  Now print out details of non-linear parts of model structure and plot
%  the shape of these non-linearities.fignu

% First check to see if there are any input non-linearities to urdisplay
flag = 0; for k=1:nu flag=flag+~strcmp(lower(G.in(k).type),'linear'); end;
plotnum=1;  % Initialise count of figure, and plot within figure

if flag
 %fignum=figure;
 %fignum=figure('Visible', 'off', 'IntegerHandle','off');
 if nargout == 1;
  fignum=figure('Visible', 'off', 'IntegerHandle','off');
  handle = [handle fignum];
 else
  fignum=figure;
 end
 %keyboard
 urdisp(['----------------------------------------------------------'],gui,guih,popup)
 urdisp(['Input Non-linearity Parameters and standard deviations:   '],gui,guih,popup)
 urdisp(['----------------------------------------------------------'],gui,guih,popup)
 for k=1:nu  % Loop over all inputs
  if ~strcmpi(G.in(k).type,'linear')
   eval(sprintf('subplot(%d1%d)',flag,plotnum));plotnum=plotnum+1;
   hold on;
   app='and standard dev:';
   urdisp([sprintf('Input block #%d of type ',k),G.in(k).type,' has estimates ',app],gui,guih,popup);
   if (strcmpi(G.in(k).type,'saturation') | strcmpi(G.in(k).type,'deadzone'))
    urdisp([' '],gui,guih,popup);
    app1 = [',sd = ',sprintf('%5.4f',SD.in(k).upper)];
    app2 = [',sd = ',sprintf('%5.4f',SD.in(k).lower)];
    urdisp(['upper limit = ',sprintf('%5.4f',G.in(k).upper),app1],gui,guih,popup)
    urdisp(['lower limit = ',sprintf('%5.4f',G.in(k).lower),app2],gui,guih,popup)
    urdisp([' '],gui,guih,popup);
    if G.in(k).upper>0 maxu = 2*G.in(k).upper; else maxu = 0.5*G.in(k).upper; end;
    if G.in(k).lower>0 minu = 0.5*G.in(k).lower; else minu = 2*G.in(k).lower; end;
    utest = minu:(maxu-minu)/100:maxu;
    if strcmpi(G.in(k).type,'saturation')
     X = sat(utest,G.in(k).lower,G.in(k).upper,1);
    else
     X = dzone(utest,G.in(k).lower,G.in(k).upper);
    end;
    h=plot(utest,X);  set(h,'Linewidth',2);
    grid; title(sprintf('Input #%d Nonlinearity',k));
    xlabel('input u'); ylabel('X(u)');
   elseif strcmpi(G.in(k).type,'hinge')
    urdisp(['eta = ' num2str(G.in(k).eta(:)','%0.4f ')],gui,guih,popup)
    if G.in(k).eta(1)>0,
     urdisp(['SD  = ' num2str(SD.in(k).eta(:)','%0.4f ')],gui,guih,popup)
    else
     urdisp(['SD  =  ' num2str(SD.in(k).eta(:)','%0.4f  ')],gui,guih,popup)
    end
    urdisp([' '],gui,guih,popup);
    if all(abs(G.in(k).eta)>0) % Avoid divide by zero when guessing relevant ranges
     breakpoints = G.in(k).eta(2:2:end)./G.in(k).eta(1:2:end);
     minu = -2*min(breakpoints);
     maxu = -2*max(breakpoints);
    else
     minu=-10; 
     maxu=10; 
    end;
    
    
    utest = minu:(maxu-minu)/100:maxu; X = hinge(utest,G.in(k).eta);
    h = plot(utest,X); set(h,'Linewidth',2);
    grid; title(sprintf('Input #%d Nonlinearity',k));
    xlabel('input u'); ylabel('X(u)');
    eta = [G.in(k).eta(:),SD.in(k).eta(:)]';
   elseif strcmpi(G.in(k).type,'poly')
    utest = -3:0.01:3; X = polynom(utest,G.in(k).eta);
    h = plot(utest,X); set(h,'Linewidth',2);
    grid; title(sprintf('Input #%d Nonlinearity',k));
    xlabel('input u'); ylabel('X(u)');
    eta = G.in(k).eta(:)';
   end;
  end; % End of check on whether we hit a linear block as we loop over all inputs
 end; % End of loop over inputs
 fignum = fignum+1;   % Subsequent plots on a new figure please
end;  % Test on whether there are any input non-linearities to urdisplay

if ~strcmpi(G.out.type,'linear')
 if nargout == 1;
  fignum=figure('Visible', 'off', 'IntegerHandle','off');
  handle = [handle fignum];
 else
  fignum=figure;
 end
 urdisp(['----------------------------------------------------------'],gui,guih,popup)
 urdisp(['Output Non-linearity Parameters and standard deviations:   '],gui,guih,popup)
 urdisp(['----------------------------------------------------------'],gui,guih,popup)
 urdisp([sprintf('Output block of type ',k),G.out.type,' has estimates and standard dev:'],gui,guih,popup);
 if (strcmpi(G.out.type,'saturation') || strcmpi(G.out.type,'deadzone'))
  urdisp([' '],gui,guih,popup);
  urdisp(['upper limit = ',sprintf('%5.4f',G.out.upper),', sd = ',sprintf('%5.4f',SD.out.upper)],gui,guih,popup)
  urdisp(['lower limit = ',sprintf('%5.4f',G.out.lower),', sd = ',sprintf('%5.4f',SD.out.lower)],gui,guih,popup)
  urdisp([' '],gui,guih,popup);
  if G.out.upper>0 maxu = 2*G.out.upper; else maxu = 0.5*G.out.upper; end;
  if G.out.lower>0 minu = 0.5*G.out.lower; else minu = 2*G.out.lower; end;
  utest = minu:(maxu-minu)/100:maxu;
  if strcmpi(G.out.type,'saturation')
   X = sat(utest,G.out.lower,G.out.upper,1);
  else
   X = dzone(utest,G.out.lower,G.out.upper);
  end;
  h = plot(utest,X);  set(h,'Linewidth',2);
  grid; title('Output Nonlinearity');
  xlabel('output y'); ylabel('X(y)');
 elseif strcmpi(G.out.type,'hinge')
  if all(abs(G.out.eta)>0) % Avoid divide by zero when guessing relevant ranges
   minu = -2*G.out.eta(3)/G.out.eta(4);
   maxu = -2*G.out.eta(length(G.out.eta)-1)/G.out.eta(length(G.out.eta));
  else minu=-10; maxu=10; end;
  
  
  minu = minu-abs(minu);
  maxu = maxu+abs(maxu);
  
  if abs(minu-maxu)>sqrt(eps),
   utest = minu:(maxu-minu)/100:maxu;
  else
   utest = linspace(-10,10,1000);
  end
  X = hinge(utest,G.out.eta);
  h = plot(utest,X); set(h,'Linewidth',2);
  grid; title('Output Nonlinearity'); xlabel('output y'); ylabel('X(y)');
  eta = [G.out.eta(:),SD.out.eta(:)]';
 elseif strcmpi(G.out.type,'poly')
  utest = -3:0.01:3; X = polynom(utest,G.out.eta);
  h = plot(utest,X); set(h,'Linewidth',2);
  grid; title('Output Nonlinearity'); xlabel('ouput y'); ylabel('X(y)');
  eta = G.out.eta(:)';
 end;
end;  % Test on output type being linear

urdisp(urline,gui,guih,popup);
urdisp(' ',gui,guih,popup);

function urdisp(str,guiRunning,guihandle,popup)

if guiRunning
 popup.append(str);
else
 disp(str);
end


