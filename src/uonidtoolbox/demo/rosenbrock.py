# function [cost,pe,grad,phi] = rosenbrock(Z,x,OPT,M,div);

# %Compute cost

# pe    = [10*(x(2)-x(1)^2) ; 1-x(1)];
# cost  = pe'*pe/2;

# if div,
#     phi  = [-20*x(1) 10 ; -1 0];
#     grad = phi'*pe;
# end

def rosenbrock(x):
    # pe = np.ndarray([10.0*(x)])

    cost = 0
    return cost