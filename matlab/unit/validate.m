%  This function is used to compare the predictive performance of an
%  estimated model G on an observed data set Z.  This data set may, or
%  may not, be the one that was used to produce G.
%
%  Usage is
%
%  [handle,cost,pe] = validate(Z,G1,G2,...);
%
%  where
%
%  Z          = [y(:),u(:)] = observed input output data
%  G1,G2,..   = Data structures which specify estimated models.
%               There may be more than one in case user wants to 
%               compare models.
%  handle     = Handle to a graphics object containing figures
%  cost       = value of quadratic cost V_N(M)
%  pe         = vector showing difference between observed output and
%               k-step ahead predictor based upon G.
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                     		  Australia.

% Copyright (C) Brett Ninness.

function [handle, cost, pe] = validate(varargin)

%Make sure we have enough input arguments
lg=length(varargin);
if lg<2, error('Need at least data Z and model M!'); end

Z=varargin{1};

%Make sure Z is in proper format
Z=startZ(Z);

%Set empty variable that will contain handles to figures
handle = [];

% Extract input and output data
[y,u,ny,nu,N] = Z2data(Z);

if nargout<2
 %Set default colour and linestyle order
 col=['b','r','k','g','m','c'];  % Set default set of colours
 col1=['b','b','k','g','m','c'];  % Set default set of colours 
 lin={'-','-.','--',':'};
 lin1={'-','-.','--',':'}; 
end;

lgd{1} = 'Data';  % Begin building legend for Plots

for i=2:lg % Loop over all models being validated

 M=varargin{i};  % Current model being validated
 
 % Bug out if model is non-parametric, we can't validate that
 if strcmpi(M.type,'nonpar'),
  handle = []; cost   = 0; pe     = 0;
  return;
 end
  
 % Get prediction error residuals
 try, 
  [cost,pe] = feval(M.costfcn,Z,m2theta(M),startOPT([]),M,0);
 catch  
  [cost,pe] = feval(M.costfcn,Z,M.theta,startOPT([]),M,0);
 end
 if size(y)~=size(pe), pe=pe'; end
 yp = y-pe;  Nyp = size(yp,1);
 
 % Initialise strings for legends on plots
 lgd{i} = ['Model ',sprintf('%d',i-1)];    
 lgda{i-1} = ['Resid. Corr. Model ',sprintf('%d',i-1)];
 lgdb{i-1} = ['Rue Model ',sprintf('%d',i-1)]; 
 if i==nargin 
  lgda{i} = 'Model 1 99% Whiteness Confidence'; 
  lgdb{i} = 'Model 1 99% Whiteness Confidence';  
 end; 
 
 if nargout<2  % Only display results if there are no output arguments requested.

  % Get colour and linestyle for plotting results for current model
	 c = col(mod(i-1,length(col))+1); l = lin{mod(i-1,length(lin))+1}; 
	 c1 = col1(mod(i-1,length(col))+1); l1 = lin1{mod(i-1,length(lin))+1}; 
  
  if i<3  % Only get new figure information for the first model, then resuse after that
   if nargout == 1; handle = [handle figure('Visible', 'off', 'IntegerHandle','off')]; else figure; end;
   startfig=gcf;
  end; 
  
  %  Display observations and 1 step ahead predictions of estimated model 
  for k=1:ny
   if (i>2) figure(startfig); end;  
   subplot(ny,1,k);
   if (i<3) plot(Z.t(1:Nyp),y(:,k),'-b','linewidth',2); end; 
   hold on;
   plot(Z.t(1:Nyp),yp(:,k),[c,l],'linewidth',2);
   legend(lgd); 
   grid;
   str = ['Observed data vs predictions of model: Output ',int2str(k)];
   title(str,'Fontsize',14);
   xlabel('Time (s)','Fontsize',14);
   set(gca,'Fontsize',14);  
  end;
  
  %  Show Correlation of residuals
  if i<3  % Only get new figure information for the first model, then resuse after that
   if nargout == 1; handle = [handle figure('Visible', 'off', 'IntegerHandle','off')]; else figure; end;
  else
   figure(startfig+1);
  end;
    
  for k=1:ny
   subplot(ny,1,k);
   % Compute sample covariance of prediction error residual on k'th output
   maxlag = min(length(pe(:,k)),26);
   R = real(ifft(abs(fft(pe(:,k)).^2)))/length(pe(:,k)); R = R(1:maxlag);
   k1 = 0:1:length(R)-1; k2 = 1:1:length(R)-1;  wun = ones(size(k2));
   if (i==2)  % Show whiteness bounds only for first model
    whitebound = 2.58/sqrt(N);  % 2.58=>99%,1.96 => 95% confidence via standard Normal
   end;
   hold on;
   plot(k1,R/R(1),c1,'linewidth',2);    
   grid;
   if (i==nargin)
    plot(k2,whitebound*wun,'-.r',k2,-whitebound*wun,'-.r','linewidth',2);
   end;   
   legend(lgda);  
   str = ['Residual Correlation + Confidence Interval: Output ',int2str(k)];
   title(str,'Fontsize',14);
   xlabel('Lag','Fontsize',14);
   set(gca,'Fontsize',14);     
  end;
  
  if nu>0 %  Show Cross-Correlation between residuals and input, but only if there is an input
    
   idx=1;
   if i<3  % Only get new figure information for the first model, then resuse after that
    if nargout == 1; handle = [handle figure('Visible', 'off', 'IntegerHandle','off')]; else figure; end;
   else
    figure(startfig+1+idx);
   end;
       
   for m=1:nu
    for k=1:ny
     subplot(ny,nu,idx); idx=idx+1;
     % Compute sample cross-covariance between pred error k and input m
     Ruep = real(ifft(fft(pe(:,k)).*conj(fft(u(:,m)))))/length(pe(:,k));
     Ruen = real(ifft(fft(u(:,m)).*conj(fft(pe(:,k)))))/length(pe(:,k));
     Rue = [flipud(Ruen(1:maxlag));Ruep(2:maxlag)]; Rue=Rue(:)';
     % Compute auto-covariance of input m and pred error k
     Ru = real(ifft(abs(fft(u(:,m)).^2)))/length(u(:,m));
     Re = real(ifft(abs(fft(pe(:,k)).^2)))/length(pe(:,k));
     sue = Ru(1)*Re(1) + 2*Ru(:)'*Re(:);  % sum_|k|<M Ru(k)Re(k)
     rue = Rue./sqrt(Ru(1)*Re(1));
     k1 = -maxlag+1:1:maxlag-1; wun = ones(size(k1));
     if (i==2)  % Figure out whiteness interval for only the first model
      whitebound = sqrt(2.58*sue/(Ru(1)*Re(1)*N));  % 2.58=>99%,1.96 => 95% confidence via standard Normal
     end; 
     hold on
     stem(k1,rue,c1,'linewidth',2);
     if (i==nargin)
      plot(k1,whitebound*wun,'-.r',k1,-whitebound*wun,'-.r','linewidth',2);
     end;
     legend(lgdb); 
     grid;
     str = ['Cross Cov Rue + CI: Output ',int2str(k),', Input ',int2str(m)];
     title(str,'Fontsize',14);
     xlabel('Lag','Fontsize',14);
     set(gca,'Fontsize',14);      
    end;
   end;
  end; % End of test on nu>0
 end; % Test on whether there are any output arguments requested
end; % Loop over number of models

hold off  % Just to make sure 