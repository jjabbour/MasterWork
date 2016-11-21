% compute the slope of the secant line of two fixed point
function [t_u,t_w] = CM_secant_line(X,k,Nc,K)
% Input:
% X = [u;c;nu]  N+2 x p
% u             N x 1
% c             1 x 1
% nu            1 x 1 
if nargin < 4
    K = 32;
    if nargin < 3
        Nc = 24;
    end
end

DIM = size(X);
N = DIM(1);
if N == (4*(Nc+1)*(2*K+1) + 1)
    t_u = X(1:N-1,k)-X(1:N-1,k-1);

elseif N == (4*(Nc+1)*(2*K+1) + 2)
    t_u = X(1:N-2,k)-X(1:N-2,k-1);
    t_w = X(end-1,k)-X(end-1,k-1);
    
end
end
