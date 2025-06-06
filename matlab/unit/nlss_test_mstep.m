function G = nlss_test_mstep(Z,M,O)

%call argmin with the correct arguments
[theta,cost,G] = argmin(Z,'VN_nlss_test',M.theta,O,M);

G.theta = theta;
G.cost  = cost;
