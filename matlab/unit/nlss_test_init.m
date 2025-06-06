function G = nlss_test_init(Z,M,O);

%Make sure to copy the input
G       = M;

% Truth is
%
% theta = [0.5; 25; 8; 0.05; sqrt(0.1); sqrt(0.1)];

G.theta = [0.4; 29; 7; 0.0; sqrt(0.1); sqrt(0.1)];
%G.theta = [0; 0; 0; 0.0; sqrt(1000); sqrt(0.1)];