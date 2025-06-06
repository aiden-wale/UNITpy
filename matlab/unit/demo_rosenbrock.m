% This test file is intended to test argmin against a well-known
% optimisation problem, i.e. Rosenbrock's banana function

clear; close all;
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end

if dsp,
	Z.y = [0;0]; % i.e. no data
	Z.u = [0;0]; % i.e. no data
	M   = []; % i.e. no model structure
	OPT = []; % i.e. use default opt settings

	n   = 5; %Number of initial points
	%x0  = repmat([1;1],1,n) + 2*randn(2,n); % initial guess
	x0  = 3.0*[...
		kron(linspace(-1,1,n),ones(1,n));...
		kron(ones(1,n),linspace(-1,1,n))];

	OPT.dsp=dsp;
	OPT.saveit=1;
	OPT.miter=100;
	%OPT.subtol=1;

	disp(' ')
	disp(' Rosenbrocks Banana Function ')
	disp(' The following optimization methods are available ')
	disp(' ')
	disp('(1) Robust Gauss-Newton')
	disp('(2) Levenberg-Marquardt')
	disp('(3) Trust Region')
	disp('(4) Quasi-Newton + trust region (QNTR)')
	disp('(5) Quasi-Newton + line search  (BFGS)')
	disp('(6) Steepest descent')
	disp(' ')
	inp=input('PLEASE SELECT ONE:  ','s');

	switch inp,
		case '1'
			OPT.dir='rgn'; OPT.subtol=1;
		case '2'
			OPT.dir='lm';
		case '3'
			OPT.dir='trust';
		case '4'
			OPT.dir='bfgs_trust';
		case '5'
			OPT.dir='bfgs';
		case '6'
			OPT.dir='grad';
		otherwise
			error('Not a valid selection');
	end


	for i=1:n*n,
		[x(:,i),cost_log,M] = argmin(Z,'rosenbrock',x0(:,i),OPT,M); Ms{i}=M;
	end

	if dsp,
		x0
		x
	end

	%Plot some stuff
	x1max=max(x0(1,:)); x1min=min(x0(1,:));
	x2max=max(x0(2,:)); x2min=min(x0(2,:));

	np=100;
	x1=linspace(min(-1.5,x1min),max(2.5,x1max),np);
	x2=linspace(min(-1.5,x2min),max(2.5,x2max),np);

	for i=1:np,
		for j=1:np,
			J(i,j)=rosenbrock([],[x1(i);x2(j)],[],[],0);
		end
	end

	if dsp,
		contour(x2,x1,log10(J));
		hold on;

		for i=1:n*n,
			plot(Ms{i}.thetait(2,:),Ms{i}.thetait(1,:),'-k*');
		end
		
		plot(1,1,'xr','linewidth',2)
	end

	echo off;
end

if dm
	disp('  ')
	disp('---------------------------------------------------------------------')
	disp('  ')
	disp('You now have access to the MATLAB workspace so that you may examine')
	disp('the results of this simulation.  To return to the demos, type "dbcont"')
	disp(' ')
	keyboard;
end;
