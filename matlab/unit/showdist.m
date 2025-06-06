%  SHOWDIST: Function to plot marginal distributions of parameters.
%  It is intended for use with the associated function POSTDIST,
%  which computes the marginals that this function plots.
%
%   Usage Example:     
%
%               z.y=y; z.u=u; m.A=4;  % Specify data and model structure
%               g=est(z,m);           % Estimate a 4th order model
%               p=postdist(z,g);      % Compute posterior dist of parameters
%               showdist(p);          % Display them
%  
%
%   Written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%             		              Australia.

%   Copyright (C) Brett Ninness

function showdist(varargin)
 
% Check number of systems to have error bounds displayed
lg=length(varargin);
if lg<1, error('Need at least one input argument'); end

col=['b','k','r','g','m','c'];  % Set default colour order
lin={'-','-.','--',':'};        % Set default linestyle order

fignum = gcf;          % Don't overwrite currently open figures
origf  = fignum;      % Remember original figure number 

for i=1:nargin,   % Loop over all systems entered
 g=varargin{i};
 nu = size(g.A,1);           % Figure out the number of inputs 
 th  = m2theta(g);           % Get standard deviations of parameter estimates
 [dummy,sd] = theta2m(th,g);
 c = col(mod(i-1,length(col))+1);  % Colour for this system plot
	l = lin{mod(i-1,length(lin))+1};  % Linestyle for this system plot
 for r=1:nu  % Look at all inputs
  fignum = origf;      
  for k=1:g.nA   % First do densities on denominator co-efficients
   %subplot([num2str(g.nA),'1',num2str(k)])
   subplot(g.nA,1,k)
   if isfield(g,'TH')   % Has the model been obtained by an MCMC run?
    plot(g.pa(k).x,g.pa(k).p,[c,l],'linewidth',2);    
    title(['Marginal posterior distribution for a_',num2str(k)]);
   else
    s=sd.A(k); v=s^2; 
    a = g.A(r,k+1); x = a-4*s:8*s/200:a+4*s; 
    N = exp(-(0.5/v)*(x-a).^2)/sqrt(2*pi*v);
    hold on; plot(x,N,[c,l],'linewidth',2); hold off;    
   end;  % Test on existance of g.TH
   ylabel('probability');
   xlabel(['a_',num2str(k)]);
   set(gcf,'Name',['Densities for Denominator Parameters for input',num2str(r)])   
   grid on;
  end;  % Loop over co-efficients in G.A
  %h=figure(fignum); fignum=fignum+1;
  
  figure
  
  for k=1:g.nB+1  % Now do densities on numerator co-efficients
   %subplot([num2str(g.nB+1),'1',num2str(k)])
   subplot(g.nB+1,1,k);   
   if isfield(g,'TH')   % Has the model been obtained by an MCMC run?
    plot(g.pb(k).x,g.pb(k).p,[c,l],'linewidth',2);    
    title(['Marginal posterior distribution for b_',num2str(k-1)])
   else
    s=sd.B(k)^2; v=s^2; 
    b = g.B(r,k); x = b-4*s:8*s/200:b+4*s; 
    N = exp(-(0.5/v)*(x-b).^2)/sqrt(2*pi*v);
    hold on; plot(x,N,[c,l],'linewidth',2); hold off;    
   end;
   ylabel('probability');
   xlabel(['b_',num2str(k-1)]);
   set(gcf,'Name',['Densities for Numerator Parameters for input',num2str(r)])   
   grid on;
  end;  % Loop over co-efficients in G.B
 end;  % Loop over number of inputs
end; % Loop over systems entered
 

