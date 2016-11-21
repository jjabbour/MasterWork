% --------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------
% September 28 2005 by Peichun Amy Tsai, Dept. of Physics, University of Toronto
%%% To calculate an appropriate exponentially decay constant, k_val, of 3Delectric potential ~ exp(-kz)
% To find out the nonlocal relation between \hat{q}_m (r) and \hat{\psi_2m} (r). 
% Last update: Oct. 17 2005
% --------------------------------------------------------------------------------------------------
% The function cheb(Nc) in this program is written by Lloyd N. Trefethen.
% --------------------------------------------------------------------------------------------------
function [Tq_vs_psi2] = qm_vs_Psi2m(K, Nc, alpha)

% Input:  K = the number of the Fourier modes from -K, .., 0, .. , K.
%             Nc = the number of the Chenbyshev collocation grids.
%             alpha = the radius ratio.
% Output: [Tq_vs_psi2] is the nonlocal relation between \hat{q}_m (r_i) and \hat{\psi_2m} (r_i).

[D,DD,x] = ccheb(Nc); % calculating the differentiation matrix for Chebyshev collocation unkowns; x = Chebyshev grids
dx = - diff(x);
delta_x = [dx(1)/2;(dx(1:end-1)+dx(2:end))/2; dx(end)/2];

Nk = 50000;
k_decay_value = 50; % k =50 the best; the decay constant k of 3d electric potential in z direction.
dk = k_decay_value/Nk; % Peichun liked k=.001...

Tq_vs_psi2 = zeros(Nc+1, Nc+1, K+1);

r = x/2 + ( alpha/(1-alpha)+1/2 ); % this is r, with r(1) = r_out and r(Nc+1) = r_in
dr = - diff(r); % this is dr
% these are the weights for an integral with respect to r.
delta_r = [dr(1)/2; (dr(1:end-1)+dr(2:end))/2; dr(end)/2];

k = [0:dk:k_decay_value]'; % k = [0:dk: k_val()]';
dk = diff(k);
delta_k = [dk(1)/2; (dk(1:Nk-1)+dk(2:Nk))/2 ; dk(Nk)/2];

[rr, kk] = meshgrid(r, k); % Matrix (K, rho') for calculating A_kr = sum( dr' r' J_m(k*r') )
[kk0, rr0] = meshgrid(k, r); % Matrix (rho, K) for calculating q_sumk_rk 
[delta_kk0,rr0] = meshgrid(delta_k,r);
[rdx, kk] = meshgrid(r.*delta_x, k); % Matrix (K, rho')

for m = 0:1:K % for m = 0:1:K, for each m
    Tq_vs_psi2(:,:,m+1) = (delta_kk0.*kk0.^2.*besselj(m,kk0.*rr0))*(rdx.*besselj(m,kk.*rr));
end

