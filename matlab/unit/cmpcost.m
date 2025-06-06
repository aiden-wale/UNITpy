function c = cmpcost(ejw,de)

N = length(ejw);
n = length(de);
v = ejw(ones(n,1),:)-de(:,ones(1,N));
c = -sum(log(conj(v(:)).*v(:)));