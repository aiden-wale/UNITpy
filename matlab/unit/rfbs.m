%  RFBS:  Solve AX=B for x when A is either lower or upper triangular by
%  using forward substitution or backward substition.
%
%  Usage is 
%
%   x=rfbs(A,B,low)
%
%  Where
%
%  A,B,X  = are the elements of the linear equation set AX=B
%
%    low  = 1 for A lower triangule, 0 fro A upper triangular
%
%   written by Brett Ninness, School of EE & CS
%              Adrian Wills   University of Newcastle
%        		              Australia.

% Copyright (C) Brett Ninness

function X = rfbs(B,A,uplow)

%uplow = 1 for lower triangular, 0 for upper

deps=100*eps; [n,m]=size(B); X=zeros(n,m);
if uplow     % Lower triangular A => Forward substition
 for i=1:m,
  if abs(A(1,1))>deps,
   X(1,i)=B(1,i)/A(1,1);
  end
  for j=2:n,
   if abs(A(j,j))>deps,
    X(j,i)=(B(j,i)-A(j,1:j-1)*X(1:j-1,i))/A(j,j);
   end
  end
 end
else         % Upper triangular A => Backward substition
 for i=1:m,
  if abs(A(n,n))>deps,
   X(n,i)=B(n,i)/A(n,n);
  end
  for j=n-1:-1:1,
   if abs(A(j,j))>deps,
    X(j,i)=(B(j,i)-A(j,j+1:n)*X(j+1:n,i))/A(j,j);
   end
  end
 end
end