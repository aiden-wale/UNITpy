%   Function to compute the probability updates
%
%   x_{t+1} ~ p_1(x_{t+1}|x_t)
%       y_t ~ p_2(y_t|x_t)
%
%   for the linear time invariant model
%
%    x_{t+1} = Ax_t + Bu_(t-delay) + Lw_t
%    y_t = Cx_t + Du_(t-delay) + v_t
%
%   Where w_t and v_t are zero mean vector white Gaussian processes with covariances
%
%   E(w_tw_t^T) = Q, E(v_tv_t^T) = R, E(w_tv_t^T) = S
%
%   and y_t, u_t are observed output and input processes.  This function
%   is intended for possible use with the SIR.m sequential importance
%   resampler routine.
%
%
%   Usage is:  [q,yp,X2] = ssmods(mode,Z,M,OPT,X,X1);
%
%   where
%   
%
%     mode==0:  the function returns the state dimension in 
%               q (i.e. q = state dimension), all other output
%               arguments are ignored.
%
%     mode==1:  the returned value q is an array of probabilities:
%
%               q(i) = p_2(y_t|x_t^i)
%
%               where x_t^i is specified by the i'th column of X.  
%               That is, it is the i'th particle at time t.
%
%               yp is a column vector of expected model output(s) for 
%               the given particles X, i.e.
%
%               yp = E{ y_t | X }
%
%               output argumnet X2 is ignored.
%
%     mode==2:  this should have exactly the same functionality as for
%               mode 1, with the addition that the returned matrix X2 
%               are sample realisations of particles x_{t+1}^i according 
%               to the model 
%
%               x_{t+1} ~ p_1(x_{t+1}|x_t)
%
%               That is, for each column of X, and with the i'th 
%               one representing a realisation of the particle x_t^i, 
%               a new particle x_{t+1}^i is formed by drawing from the
%               probability distribution 
%
%               x_{t+1}^i ~ p_1(.|x_t^i).
%
%     mode==3:  the returned value q is an array of
%               probabilities:
%
%               q(i) = p_1( x^i_{t+1|N} | x^i_{t|t} )
%
%               where x^i_{t+1|N} is the i'th column of X1 and
%               x^i_{t|t} is the i'th colummn of X.
%
%               yp(:,i) should contain the prediction x^i_{t+1|t}
%
%               output argument X2 is ignored.
%
%     mode==4:  the returned variable q is an array of
%               probabilities:
%
%               q(i) = p_1( x^i_{t+1|N} | x^i_{t+1|t} )
%
%               where x^i_{t+1|N} is the i'th column of X1 and
%               x^i_{t+1|t} is the i'th colummn of X (provided 
%               by the previous call with mode==3)
%
%               yp and X2 output arguments are ignored.
%
%   Z:          Input-Output data in one of two forms.  The standard form 
%               is for it to be a record with elements Z.y and Z.u, each
%               of which are matrices with number of rows equal to the
%               number of data samples, and number of columns equal (respectively)
%               to the number of outputs and the number of inputs.  On
%               the other hand, Z can be a matrix of the form Z = [y,u] 
%               where it is assumed that y is a column vector of output 
%               measurements and u is a matrix whose columns are the
%               input measurements; in this latter MISO models are 
%               being considered.  
%
%   M:          Data structure which defines the above model:
%
%               As an example of how to write func, see ssmods.m which
%               implements the LTI model detailed above.
%
%
%   OPT:        Data structure which defines options for the estimation
%               algorithm as follows:
%    OPT.pnum   Number of particles to use.  Default is OPT.pnum=100.
%
%
%                             

% Copyright (C) Brett Ninness

function [varargout]=ssmods(mode,Z,M,OPT,X,X1)
 
nx = size(M.ss.A,1);  % Determine state dimension

if mode==0, 
    varargout(1)={nx};
    return;
end

if mode==1 | mode==2,
    % Compute Measurement update weights q_t^i = p(y_t|x^i_t|t-i)
    yp = M.ss.C*X(:,:)+M.ss.D*repmat(Z.u(:),1,OPT.pnum);
    pe=repmat(Z.y(:),1,OPT.pnum)-yp;
 
 PY=pe*pe'/OPT.pnum;
 
    qpe = pe.*(M.ss.R\pe);
    if size(qpe,1)>1 qpe=sum(qpe); end;
    q=exp(-qpe/2); q=q/sum(q);        

    % Pass back predicted y based on passed in (predicted) state
    yp=mean(yp,2);

    varargout(1)={q}; varargout(2)={yp};
 
    varargout(3)={PY};
 
    if mode==2,
        % Now do time update of particles x^i_t+1|t \sim p(x_t+1|x_t^i)
        X=M.ss.A*X+M.ss.B*repmat(Z.u(:),1,OPT.pnum)+sqrtm(M.ss.Q)*randn(nx,OPT.pnum);
        varargout(3)={X};
    end
    return;
end

if mode==3,
    % Perform prediction x_{t+1|t}
    X=M.ss.A*X+M.ss.B*repmat(Z.u(:),1,OPT.pnum);
    pe=X1-X; qpe=sum(pe.*(M.ss.Q\pe),1);
    q=exp(-qpe/2); q=q/sum(q);
    varargout(1)={q}; varargout(2)={X};
end

if mode==4,
    % Calc p_1(x_{t+1|N} | x^i_{t+1|t})
    pe=X1-X; qpe=sum(pe.*(M.ss.Q\pe),1);
    q=exp(-qpe/2);
    q=q/sum(q);
    varargout(1)={q};
end
