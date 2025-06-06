function G=em_hamm(Z,M,OPT);

% Extract sizes of input and output from data matrix
[y,u,ny,nu,N] = Z2data(Z);

G=M;            % Initial estimate becomes candidate for final estimate
nx=size(G.ss.A,1); % Get state dimension
stop_crit = 0;  % Flag for termination of search
count     = 0;  % Count on number of search iterations

if OPT.dsp,     % Let people know the algorithm being used.
	disp('EM based Hammerstein System Estimation');
	disp('--------------------------------------');
end;

% Initialize uncertainty in Hammerstein system components
mz.ss.Q  = 10*eye(nu,nu);
mz.ss.P1 = eye(nu,nu);
mz.ss.X1 = zeros(nu,1);

G.ss.X1=zeros(nx,1);
G.ss.P1=eye(nx);
G.ss.Q=100*eye(nx);
G.ss.R=0.01*eye(ny);
G.ss.S=zeros(nx,ny);

% Now loop over E and M steps
while stop_crit==0,
	count = count + 1;

	% Perform the E-step using a Kalman smoother
	m=G; 
	opt.allP=1; 
	m.in=G.in;
	z=u2x(u,m);  % Pass input through existing estimate of input non-linearity

	% Form Augmented state space model that includes Hamm output as a state
	mm.ss.A=[m.ss.A,m.ss.B;zeros(nu,nx),zeros(nu,nu)]; mm.ss.B=[zeros(nx,nu);eye(nu,nu)];
	mm.ss.C=[m.ss.C,m.ss.D]; mm.ss.D=zeros(ny,nu);
	mm.ss.X1=[m.ss.X1;mz.ss.X1]; mm.ss.P1=[m.ss.P1,zeros(nx,nu);zeros(nu,nx),mz.ss.P1];
	mm.ss.Q=[m.ss.Q,zeros(nx,nu);zeros(nu,nx),mz.ss.Q]; mm.ss.R=m.ss.R;
	mm.ss.S=[m.ss.S;zeros(nu,ny)];

	% Find smoothed estimate of system state for current parameters
	zz.y = y; zz.u = [z(2:end,:);zeros(1,nu)];
	g = ks(zz,mm,opt);
	zhat = g.ss.X(nx+1:nx+nu,1:N);  % Estimate of Hamm Non-lin output

	G.PE(count) = ny*g.mse; G.LL(count) = g.LL;
	if OPT.dsp,
		if G.LL(count) < 0,
			disp(sprintf( 'Iteration # = %3d, LL = %3.7e, PE cost = %3.7e',count,G.LL(count),G.PE(count)));
		else
			disp(sprintf( 'Iteration # = %3d, LL =  %3.7e, PE cost = %3.7e',count,G.LL(count),G.PE(count)));
		end
	end;

	% Use Results of KS step to find necessary conditional expectations
	Ptsum=zeros(nx+nu,nx+nu); Mtsum = Ptsum;
	for t=1:N Ptsum=Ptsum+g.ss.Pt{t};  Mtsum=Mtsum+g.ss.Mt{t}; end;
	xy = [g.ss.X(1:nx,2:end); y'];  Phi = xy*xy';
	Psi = xy*g.ss.X(:,1:N)';  Sigma=g.ss.X(:,1:N)*g.ss.X(:,1:N)';
	Sigma = Sigma+Ptsum; Psi(1:nx,:)=Psi(1:nx,:)+Mtsum(1:nx,:);
	p1=g.ss.Pt{1}; pn=g.ss.Pt{N+1};
	Phi(1:nx,1:nx)=Phi(1:nx,1:nx)+Ptsum(1:nx,1:nx)-p1(1:nx,1:nx)+pn(1:nx,1:nx);

	% Now for the M-Step: Re-estimate parametrization of linear dynamics first
	H = chol([Sigma Psi'; Psi Phi]);
	Gamma = H(1:nx+nu,nx+nu+1:end)'/H(1:nx+nu,1:nx+nu)';
	Pi = H(nx+nu+1:end,nx+nu+1:end)'*H(nx+nu+1:end,nx+nu+1:end)/N;
	G.X1 = g.ss.X(1:nx,1); Px = g.ss.Pt{1}; G.P1=Px(1:nx,1:nx);

	% separate out the system matrices.
	G.A = Gamma(1:nx,1:nx); G.C = Gamma(nx+1:end,1:nx);
	if nu>0, G.B = Gamma(1:nx,nx+1:end); G.D = Gamma(nx+1:end,nx+1:end); end;
	G.Q = Pi(1:nx,1:nx); G.S = Pi(1:nx,nx+1:end); G.R = Pi(nx+1:end,nx+1:end);

	% Now re-estimate parametrization of input non-linearity
	for k=1:nu  % Do once for each input
		if ~strcmp(M.in(k).type,'linear')
			m1.A=0; m1.in = G.in(k);
			opt.dsp=0; opt.subtol=1e-18;
			g1 = est([zhat(k,:)',u(:,k)],m1,opt); G.in(k)=g1.in;
		end;
	end;
	zest     = u2x(u,G);
	zerr     = zhat'-zest;
	mz.ss.Q  = (zerr'*zerr+Ptsum(nx+1:nx+nu,nx+1:nx+nu))/N;
	mz.ss.P1 = Px(nx+1:nx+nu,nx+1:nx+nu);
	mz.ss.X1 = g.ss.X(nx+1:nx+nu,1);

	%Update stopping criterion
	if (count>=OPT.miter), stop_crit=1; end
	if count>1,
		Ldiff = G.LL(count)-G.LL(count-1); Rdiff = Ldiff/(abs(G.LL(count))+abs(G.LL(count-1)));
		if (Ldiff < 0 )
			disp('Quitting due to decrease in likelihood'); stop_crit = 1;
		elseif (Rdiff < OPT.mdec)
			disp('Quitting due to rel. gain');  stop_crit=1;
		end;
	end
end; % for loop
