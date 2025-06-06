%  Running this demos estimation of a static (memoryless) 
%  non-linearity

clear; close all;  
global dm; if isempty(dm), clear global dm; dm=0; end
global dsp; if isempty(dsp), clear global dsp; dsp=1; end
global trans;
if dsp>1, clc; echo on; end
OPT.dsp = dsp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Experiment Conditions                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N   = 1000;           % Number of samples
var = 1e-4;    % White Measurement Noise variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify Static Non-linear system to be estimated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nonlin_type = 'atan';   % Can be 'sat', 'atan', 'poly', 'dzone';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulate a data record                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = 10*rand(1,N); u=u-max(u)/2*ones(size(u));  % Input
noise = sqrt(var)*randn(size(u));              % Measurement Noise
if strcmp(nonlin_type,'sat')			
 y = sat(u,-0.3,0.3,1);
elseif strcmp(nonlin_type,'dzone')			  
 y = dzone(u,-2,2);
elseif strcmp(nonlin_type,'atan')			  
 y = atan(u);
elseif strcmp(nonlin_type,'poly')			  
 y = 0.2*u + 0.3*u.^2-0.1*u.^3+0.1*u.^8;
end;
y = y(:) + noise(:);
Z = [y(:),u(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify a Model Structure to fit to the data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M.type     = 'static';
M.in.type  = 'hinge';   % Non-linearity type - could by `poly',
                        % 'saturation' or 'deadzone' as well.

if strcmp(M.in.type,'poly'); % Initial guess at parameterisation of non-linearity.
 M.in.eta = [1,zeros(1,11)];  
elseif strcmp(M.in.type,'hinge');
 M.in.eta = [0.05,1,-0.05,-1,-0.05,1];
elseif strcmp(M.in.type,'deadzone'),
	M.in.upper=0.5;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Estimate on basis of noise corrupted data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = est(Z,M,OPT); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dsp,
	if strcmp(M.in.type,'poly')
		h1=plot(u(:),[y(:),polynom(u(:),G.in.eta)],'x');
	elseif strcmp(M.in.type,'hinge')
		h1=plot(u(:),[y(:),hinge(u(:),G.in.eta)],'x');
	elseif strcmp(M.in.type,'saturation')
		h1=plot(u(:),[y(:),sat(u(:),G.in.lower,G.in.upper,1)],'x');
	elseif strcmp(M.in.type,'deadzone')
		h1=plot(u(:),[y(:),dzone(u(:),G.in.lower,G.in.upper)],'x');
	end;
	set(h1,'Linewidth',2);
	legend({'Measured',['Estimated via ' M.in.type ' function']},'location','southeast')
	grid
end

echo off;

if dm
 disp('  ')
 disp('---------------------------------------------------------------------')
 disp('  ')
 disp('You now have access to the MATLAB workspace so that you may examine')
 disp('the results of this simulation.  To return to the demos, type "dbcont"')
 disp(' ')
 keyboard; 
end;







