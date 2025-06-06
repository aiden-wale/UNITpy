function [flag,message] = checkZM(Z,M)

flag    = 0;
message = '';

Z = startZ(Z);

if Z.nu ~= M.nu,
 flag = 1;
 message  = ['Number of inputs in data (' num2str(Z.nu) ') does not equal number of inputs in model (' num2str(M.nu) ').'];
end

if Z.ny ~= M.ny,
 flag = 1;
 message  = ['Number of outputs in data (' num2str(Z.ny) ') does not equal number of outputs in model (' num2str(M.ny) ').'];
end

