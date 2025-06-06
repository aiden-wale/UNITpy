%   EMNLSS: Subroutine called by EST that computes
%   expectation-maximisation steps for state space model 
%   structures.  This routine is not meant to be directly called by the
%   user - call est.m instead which sets up all the information that needs
%   to be passed to EMNLSS.m
%
%    written by Brett Ninness,   School of EE & CS
%               Adrian Wills,    University of Newcastle
%                                Australia.

% Copyright (C) Brett Ninness, Adrian Wills


% EM routine for NLSS models
function G = emnlss(Z,M,OPT)


% Unspecified parts of OPT -> defaults
if ~exist('OPT'), 
 OPT = startOPT([]); 
else
 OPT = startOPT(OPT); 
end

% Copy EM algorithm result into output structure G
G=M; % Input specs on model type etc are copied to output

% Call initialisation routine - focus is to initialise parameter vector G.theta
G = feval(M.nlss.init,Z,G,OPT);

% Save the thetas
thetait = zeros(length(G.theta(:)),OPT.emit+1);
thetait(:,1) = G.theta(:);

% Start Main routine
for i=1:OPT.emit,
 if OPT.dsp disp(sprintf('Iteration #%4i',i)); end;
 
 %Call the E-step - call particle smoother to set G.xs
 G = feval(M.nlss.estep,Z,G,OPT);
 
 %Call the M-step - gradient based search (or whatever you like) to
 %minimize Q - this gives new G.theta
 G = feval(M.nlss.mstep,Z,G,OPT);
 
 thetait(:,i+1) = G.theta(:);
end

% Add legend for prospective plotting
G.disp.legend=['Estimated ',upper(G.type),' model via EM'];

% Record that EM algorithm was used
G.alg='em';

% Save the theta iterations
G.thetait = thetait;