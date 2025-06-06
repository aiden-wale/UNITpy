% SYSTR - function to simulate a bilinear system.  Never meant to be
% called by the user - strictly for internal use
%
%
%   written by Brett Ninness, School of EE & CS
%                             University of Newcastle
%                             Australia.
%

% Copyright (C) Brett Ninness

function [xh,yh,pe,ukx] = systr(y,u,A,B,C,D,K,F,G,X1,isD,isK,isF,isG,isX1,mukx)

N = max(length(y),length(u));
n = size(A,1);
m = size(B,2);
p = size(C,1);

x   = zeros(1,n);
xh  = zeros(N,n);
yh  = zeros(N,p);
pe  = zeros(N,p);
ukx = zeros(N,mukx*n);
At  = A.';
Bt  = B.';
Ct  = C.';
Dt  = D.';
Kt  = K.';
Ft  = F.';
Gt  = G.';
if isX1,
    x(1,:) = X1(:)';
end
for t=1:N,
 %Save state into xh
 xh(t,:) = x;
 
 %Compute u kron x
 ukx(t,:) = kronaw(u(t,1:mukx),xh(t,:));
 
 %Update output
 yh(t,:) = xh(t,:)*Ct + u(t,:)*Dt + ukx(t,:)*Gt;
 pe(t,:) = y(t,:) - yh(t,:);
 
 %Update state
 x = xh(t,:)*At + u(t,:)*Bt + ukx(t,:)*Ft + pe(t,:)*Kt;
end