% M2SYS - function to convert an estimated model structure into Matlab
% Control Systems Toolbox sys object.
%
% Usage is:
%
% [sysG,sysH] = m2sys(G)
%
% Where
%
% G     = Estimated model structure
%
% sysG  = Input/output dynamics as a TF sys object
%
% sysH  = Noise model as a TF sys object
%
% written by Brett Ninness, School of EE & CS
%            Adrian Wills   University of Newcastle
%      		                Australia.

% Copyright (C) Brett Ninness.

function [sysG,sysH] = m2sys(g)

% Put polynomial into correct format and account for delays
sysG=[]; sysH=[];

if strcmpi(g.type,'nonpar'),
 sysG = frd(g.G,g.w);
else
 A=[]; for i=1:size(g.A,1), A{i}=g.A(i,:); end
 B=[]; for i=1:size(g.B,1), B{i}=[zeros(1,g.delay(i)) g.B(i,:)]; end
 C=g.C; D=g.D;

 switch g.type,
  
  case 'ss'
   
  if g.op=='q'
   
   %Account for any delays
   delay=g.delay;
   nx=size(g.ss.A,1);
   nu=size(g.ss.B,2);
   ny=size(g.ss.C,1);
   A=g.ss.A; B=g.ss.B; C=g.ss.C; D=g.ss.D;
   if isempty(D),
    D = zeros(ny,nu);
   end
   for i=1:nu,
    for k=1:delay(i),
     A=[A B(:,i);zeros(1,nx+1)];
     B(:,i)=zeros(nx,1);
     B=[B;[zeros(1,i-1) 1 zeros(1,nu-i)]];
     C=[C D(:,i)];
     D(:,i)=zeros(ny,1);
     nx=nx+1;
    end
   end
   sysG=ss(A,B,C,D,g.T);
   if isfield(g.ss,'K'),
    if ~isempty(g.ss.K),
     p=size(g.ss.C,1);
     sysH=ss(g.ss.A,g.ss.K,g.ss.C,eye(p),g.T);
    end
   else
    n=size(g.ss.A,1); p=size(g.ss.C,1);
    sysH=ss(zeros(n),zeros(n,p),zeros(p,n),eye(p),g.T);
   end

  elseif g.op=='s'
   
   nx=size(g.ss.A,1);
   nu=size(g.ss.B,2);
   ny=size(g.ss.C,1);
   A=g.ss.A; B=g.ss.B; C=g.ss.C; D=g.ss.D;
   if isempty(D),
    D = zeros(ny,nu);
   end
   sysG=ss(A,B,C,D);
   if isfield(g.ss,'K'),
    if ~isempty(g.ss.K),
     p=size(g.ss.C,1);
     sysH=ss(g.ss.A,g.ss.K,g.ss.C,eye(p));
    end
   else
    n=size(g.ss.A,1); p=size(g.ss.C,1);
    sysH=ss(zeros(n),zeros(n,p),zeros(p,n),eye(p));
   end
  end

  case {'fir','nfir'}
      
   if g.op=='s',
    sysG=tf(B,1,'variable','s');
    sysH=tf(1,1,'variable','s');
   else
    sysG=tf(B,1,g.T,'variable','z^-1');
    sysH=tf(1,1,g.T,'variable','z^-1');
   end

 case {'ar','nar','arma','narma'}
      
   if g.op=='s',
    sysG=[];
    sysH=tf(C,A,'variable','s');
   else
    sysG=[];
    sysH=tf(C,A,g.T,'variable','z^-1');
   end

 case {'arx','narx'}
      
   if g.op=='s',
    sysG=tf(B,A,'variable','s');
    sysH=tf([1 zeros(1,g.nA)],A{1},'variable','s');
   else
    sysG=tf(B,A,g.T,'variable','z^-1');
    sysH=tf(1,A{1},g.T,'variable','z^-1');
   end

 case {'oe','noe'}
 
   if g.op=='s',
    sysG=tf(B,A,'variable','s');
    sysH=tf(1,1,'variable','s');
   else
    sysG=tf(B,A,g.T,'variable','z^-1');
    sysH=tf(1,1,g.T,'variable','z^-1');
   end

 case {'armax','narmax'}
      
   if g.op=='s',
    sysG=tf(B,A,'variable','s');
    sysH=tf(C,A{1},'variable','s');
   else
    sysG=tf(B,A,g.T,'variable','z^-1');
    sysH=tf(C,A{1},g.T,'variable','z^-1');
   end

 case {'bj','nbj'}
      
   if g.op=='s',
    sysG=tf(B,A,'variable','s');
    sysH=tf(C,D,'variable','s');
   else
    sysG=tf(B,A,g.T,'variable','z^-1');
    sysH=tf(C,D,g.T,'variable','z^-1');
   end

 end
end