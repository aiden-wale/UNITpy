
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%                             Australia.
%

% Copyright (C) Brett Ninness

function [dy] = sysfr(N,y,XB,CX,dAKC,dBKD,dK,dC,dD,isK,isD)

if isK && isD,
 for k=1:N,
  dX            = dBKD + dAKC*XB(:,:,k) + dK*y(:,:,k);
  dy(:,:,k)     = dC*XB(:,:,k) + CX(:,:,k)*dX + dD;
 end
elseif isK && ~isD,
 for k=1:N,
  dX            = dBKD + dAKC*XB(:,:,k) + dK*y(:,:,k);
  dy(:,:,k)     = dC*XB(:,:,k) + CX(:,:,k)*dX;
 end
elseif ~isK && isD,
 for k=1:N,
  dX            = dBKD + dAKC*XB(:,:,k);
  dy(:,:,k)     = dC*XB(:,:,k) + CX(:,:,k)*dX + dD;
 end
else
 for k=1:N,
  dX            = dBKD + dAKC*XB(:,:,k);
  dy(:,:,k)     = dC*XB(:,:,k) + CX(:,:,k)*dX;
 end
end