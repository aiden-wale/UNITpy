function M = sscan_t2m(M,theta)

M.ss.A  = [theta(1) 1; theta(2) 0];
M.ss.B  = [theta(3);theta(4)];
M.ss.C  = [1 0];
M.ss.D  = [];

K = [theta(5);theta(6)];

M.ss.Q  = K*K';
M.ss.S  = K;
M.ss.R  = 1;