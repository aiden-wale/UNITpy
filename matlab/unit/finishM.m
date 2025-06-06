function G = finishM(Z,M,OPT);

G = M;  % Pass all input model structure information to output

% Get appropriate text description of model structure;
ep = estmap(Z,G,OPT); G.modelEquations = ep.modelEquations;

% Add a matlab-style system description
[G.sysG,G.sysH]=m2sys(G);

% Add a frequency response
G = m2f(G);

% Reorder fields in alphabetical order
G=orderfields(G);
