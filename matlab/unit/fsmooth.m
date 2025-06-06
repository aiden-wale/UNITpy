function [PP,Pssum,Psejw,Psej2,MSE,LL] = fsmooth(V,S,Vi,B,C,D,Q,R,w,ejw,Y,U,discrete) 


%/* Get sizes */
n  = size(V,1); %/* number of states */
m  = size(B,2);  %/* number of inputs */
mm = size(U,2);  %/* number of inputs in input data, this handles the case of FRF and FFT based data*/
p  = size(C,1);  %/* number of outputs */
N  = max(size((Y)));  %/* number of frequency points */
	
%/*---------------------------------------------*/
%/* Main computational activity starts here */
    
%/* Get Cholesky factors of Q */
qq = rchol(Q);

Q=qq;

%/* Clear some variables */
Pssum = zeros(n,n);
Psejw = zeros(n,n);
Psej2 = zeros(n,n);
PP    = zeros(2*n+m+p,N*mm);
MSE   = 0;
LL    = 0;

for k=1:N,
	%/* Construct (ejw(k)*I - A)^{-1} = U*(ejw(k)*I-S)^{-1}*V*/
	Ak  = V*diag(1./(ejw(k)-S))*Vi;

	%/* Now form Bk = A_k^{-1}B */
	Bk = Ak*(B*U(:,:,k));
	
	%/* Ak = Ak*Q^{1/2} */
	Ak = Ak*Q';
		
	%/* Ck = C*Ak = C*(ejw(k)*I-A)^{-1}Q^{1/2} */
	Ck = C*Ak;
        
	%/* Form Ek = Y - C*(ejw(k)*I-A)^{-1}*B - D */
	Ek = Y(:,:,k) - C*Bk - D*U(:,:,k);
	
	%/* Update MSE */
	MSE = MSE + Ek(:)'*Ek(:);
        
	%/* Form Rk = chol(Ck*Ck' + R) */
	Rk = chol(Ck*Ck' + R);
	

	%/* Let Ek = Rk^{-1/2}*Ek */
	Ek = Rk\Ek;

		
	%/* Update Log-Likelihood */
	alp = prod(diag(Rk));
	LL  = LL - 2*log(alp) - Ek(:)'*Ek(:);
		
	%/* Let Ck = Rk^{-'/2}*Ck */
	Ck = Rk\Ck;
		
	%/* Fk = Ak * Ck'*/
	Fk = Ak*Ck';
		
	%/* Form Xk = Bk + Fk*Ek */
	Xk = Bk + Fk*Ek;
		
	%/* Form Pk = Ak*Ak' - Fk*Fk' */
	Pk = Ak*Ak' - Fk*Fk';
		
	%/* Fill columns of PP */
	PP(:,(k-1)*mm+1:k*mm) = [ejw(k)*Xk;Y(:,:,k);Xk;U(:,:,k)];
		
	%/* Compute sums */
	Pssum = Pssum + Pk;
	Psejw = Psejw + ejw(k)*Pk;
	if ~discrete,
		Psej2 = Psej2 + w(k)^2*Pk;
	end
end

if discrete, Psej2 = Pssum; end