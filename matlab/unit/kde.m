%  KDE: Univariate probability density estimation via kernel method,
%  ie. kernel density estimation.  That is, given discrete estimates p_k
%  of the probability of event x_k (for example, obtained via a
%  histgoram), the density estimate f(x) is computed as
%
%  f(x) = 1/(n*v) \sum_{k=1}^n p_k K[ (x-x_k)/v ]
%
%  where K(x) is a smoothing kernel of unit area.
%
%  Usage is
%
%  hest = kde(h,OPT);
%
%  where:
%  h:           Data structure specifying discrete estimates p_k of
%               events x_k.  It contains elements
%  h.p          Vector of estimates p_k;
%  h.x          Vector of events x_k.  That is, h.p(k) = p_k = p(x_k).
%
%  OPT:         Data structure which defines options for the estimation
%               algorithm as follows:
%  OPT.krange:  Width of support of the kernels as an *even* multiplier of
%               the width of the bins in the input histogram.  Default is
%               OPT.krange=6.
%  OPT.knum:    Number of points per bin specified in h.x to be used
%               within OPT.krange to represent the kernels.  Default is
%               OPT.knum=6;  So, as example, with defaults, if
%               elements in h.x are spaced 0.1 units apart, then
%               OPT.krange will be 6*0.1=0.6 to span three bins either
%               side, and within this range there will be 6 points per
%               bin, or 6*6 = 36 points over the 0.6 range used to
%               specify the kernels.
%  OPT.kernel   A string specifying type of kernel K(x) to be
%               used. Options are 'gaussian' and 'triangle'. Default is
%               'gaussian';
%  OPT.v        Scaling paramter v used to expand or contract K(x) as
%               indicated above in specification of kernel density estimate.
%               Default is OPT.v = (h.x(2)-h.x(1))^2 - that is, the square
%               of the distance between "bins" given in h.x;  
%  hest.p       Vector of kernel smoothed probability denisty function
%               estimates
%  hest.x       Corresponding x axis so that plot(hest.x,hest.p) plots
%               estimated probability density function.
%  hest.pin     Input probabilities normalised so that with respect to
%               h.x, the area underneath the associated input estimate is
%               one.  This allows input and output to be compared as
%               bar(h.x,hest.pin); hold on; plot(hest.x,hest.p); hold off;
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%        		              Australia.

% Copyright (C) Brett Ninness

function hest = kde(histogram,OPT);

% For ease of reference shorten name
bin_p = histogram.p; bin_x = histogram.x;
bin_diff = bin_x(2)-bin_x(1);                          % x axis distance between bins; 

% Normalise input histogram to have unit area
bin_p=bin_p/(sum(bin_p)*bin_diff);

if ~exist('OPT') OPT.krange=6;                    end;
if ~isfield(OPT,'krange') OPT.krange=4;           end;  % Width of kernel support as integer multiple of bin_diff - must be even
if ~isfield(OPT,'knum')   OPT.dnum=4*OPT.krange;  end;  % Number of points representing kernel - must be even multipler
if ~isfield(OPT,'kernel') OPT.kernel='gaussian';  end;  % Kernel type
if ~isfield(OPT,'v')      OPT.v=(bin_diff^2); end;      % Kernel scaling (variance)

krange = OPT.krange*bin_diff;                          % Width of kernel support
del_x = krange/OPT.dnum;                               % x-axis increment between samples of estimated density 
x=-krange/2:del_x:krange/2;                            % x-axis support of kernel
del_n =  floor(bin_diff/del_x);                        % Number of x_axis samples between bins of input histogram
lenx = length(x);                                      % Just to make things more readable below

% Compute appropriate kernel
if strcmpi(OPT.kernel,'gaussian')
 f = exp(-1/2*(x.^2)/OPT.v)/sqrt(2*pi*OPT.v);
elseif strcmpi(OPT.kernel,'triangle') 
 f=max(0,(1/OPT.v-(1/OPT.v^2)*abs(x)));
else
 error('Unknown specification for OPT.kernel: type "help kde"');
end;

% Now compute smoothed estimate
hest.x = bin_x(1)-krange/2:del_x:bin_x(end)+krange/2;
hest.x = hest.x(1:floor(length(hest.x)/2)*2); % Make sure of consistency of length (Matlab bug in that sometimes last element not included?)
hest.p = zeros(size(hest.x));  
for k=1:length(bin_x)  % Loop through all bin locations
 idx = 1+(k-1)*del_n:lenx+(k-1)*del_n; 
 idx = idx(idx<=length(hest.p));  % If bin_diff is odd then idx could overflow on last frame
 hest.p(idx) = hest.p(idx) + bin_p(k)*f(1:length(idx));
end;

%Output smoothed estimate after normalisation
hest.p=hest.p/(sum(hest.p)*del_x);

% Output normalised input
hest.pin = bin_p;



