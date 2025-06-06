%  Function to look for outliers and remove them
%
%  Usage is z = outlier(x)
%
%  where x = vector containing outliers
%        z = vector with outliers identified, and replaced with median of
%            remaining values.
%      
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.

% Copyright (C) Brett Ninness, 

function z = outlier(x);

%  Look at empirical distribution of x
[a,b] = hist(x);

if(sum(a==0)>3)%  Is there a gap of more than 3 bins in the histogram?
  %  Find the bin with the outlier  
  [dummy,index]=max(abs(diff(b(a~=0))));
  
  %  get rid of outlier  
  if (index==1)  %  Outlier at left of histogram
    x(x<b(2)) = median(x);
  else
    x(x>b(length(b)-1)) = median(x);
    end;
end;

z=x(:)';  % Make sure a row vector is returned.
