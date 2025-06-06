function G = finishMctstf(Z,M,OPT);

G = M;  % Pass all input model structure information to output

% We got time domain tf model by grey box ss, but hide that detail
G.type = 'oe'; G = rmfield(G,'par'); G = rmfield(G,'t2m');
G.C = []; G.d = [];

% Get appropriate description of tf model structure
ep = estmap(Z,G,OPT); G.modelEquations = ep.modelEquations;

% Add a matlab style system description
[G.sysG,G.sysH]=m2sys(G);

% Add a frequency response
G = m2f(G);

% Compute covariance matrix of systen estimates
R = G.jacobian; R = (R'*R);  G.P = G.var*pinv(R);

% Reorder fields in alphabetical order
G=orderfields(G);