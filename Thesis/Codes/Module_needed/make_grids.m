function [D,D2,yi,I,It,theta,riprime] = make_grids(Nc,K,rr)

[D, D2,yi] = ccheb(Nc); % calculating the differentiation matrices for Chebyshev collocation unkowns; yi = Chebyshev grids
% yi is of length Nc+1 and goes from 1 to -1.
% D, and D2 are matrices of size (Nc+1)x(Nc+1) that produce 
% derivatives in the radial direction

I = eye(Nc-1, Nc-1); It = eye(Nc+1, Nc+1);

% 2K+1 grid points in the theta direction.  Includes 0, excludes 2 pi
theta = 2*pi*(0:1:2*K)/(2*K+1);       % grid points in theta direction

% there are Nc+1 points between [-1 and 1].  These are the Chebyshev grids.

%  = r'; the transformation between r = [ri/d, ro/d] and x = [-1,1];
% rr is the radius ratio (alpha in Peichun's thesis)
% rr/(1-rr) = r_in/d = r_in' (the nondimensionalized inner radius)
% 1/(1-rr) = 1 + rr/(1-rr) = r_out/d = r_out' (nondimensionalized outer radius)
%
% riprime goes from 1 + rr/(1-rr) to rr/(1-rr).  (Outer to inner)
riprime = yi/2 + ( rr/(1-rr)+ 1/2 );  


