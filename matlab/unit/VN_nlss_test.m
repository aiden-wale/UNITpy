function [cost,pe,grad,phi] = VN_nlss_test(Z,theta,OPT,ARGS,div)

%form prediction error
y    = Z.y(:)';
u    = Z.u(:)';
x    = ARGS.xs;
M    = OPT.pnum;
N    = length(y);
cost = 0;
pe   = [];
grad = [];
phi  = [];

uu   = u(ones(M,1),:);
yy   = y(ones(M,1),:);
xt   = squeeze(x(:,:,1:N));
xt1  = squeeze(x(:,:,2:N+1));
xs   = xt.*xt;
x1x2 = xt./(1.0 + xs);
e1   = (xt1 - (xt*theta(1) + x1x2*theta(2) + uu*theta(3)))/theta(5);
e2   = (yy-xs*theta(4))/theta(6);
cost = log(theta(5)^2) + log(theta(6)^2) + (sum(e1(:).*e1(:)) + sum(e2(:).*e2(:)))/M/N;

if div,
    grad    = zeros(6,1);
    tmp1    = 2.0/(M*N*theta(5));
    tmp2    = 2.0/(M*N*theta(6));
    grad(1) = -sum(xt(:).*e1(:))*tmp1;
    grad(2) = -sum(x1x2(:).*e1(:))*tmp1;
    grad(3) = -sum(uu(:).*e1(:))*tmp1;
    grad(4) = -sum(xs(:).*e2(:))*tmp1;
    grad(5) =  sum(ones(M*N,1)-e1(:).*e1(:))*tmp1;
    grad(6) =  sum(ones(M*N,1)-e2(:).*e2(:))*tmp2;    
end