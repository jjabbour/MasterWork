function [Jqphi_m,Jwphi_m,Jpsi2q_m] = compute_nonlinearity(psi2m,wm,phim,qm,Phi0_r,Q0_r,W0_r,Psi20_r,riprime,Nonlinear,D,Nc,q)

% =====================================================================
% the aliasing removal method for the pseudo-spectral method
% ================================================

% we are given spectral data at K modes.  We want to expand this to spectral data at q=2*K modes, do an IFFT, compute the nonlinearity at 4K+1 points, do a FFT, resulting in spectral data at 2K modes, and then we pass back out of the subroutine the spectral data at K modes.

K=q/2;
% Need to pad out the inputs for the dealiasing...
psi2m(:,K+2:2*K+1)=zeros(Nc+1,K);
wm(:,K+2:2*K+1)=zeros(Nc+1,K);
phim(:,K+2:2*K+1)=zeros(Nc+1,K);
qm(:,K+2:2*K+1)=zeros(Nc+1,K);
% Now that I've extended the FFT coefficients from having modes 0 ...
% K to having modes 0 ... q=2*K, I need to multiply them by a scaling
% factor because of how FFTs work...
psi2m = psi2m*(2*q+1)/(2*K+1);
wm = wm*(2*q+1)/(2*K+1);
phim = phim*(2*q+1)/(2*K+1);
qm = qm*(2*q+1)/(2*K+1);

% dphi/dtheta
testm = phim*diag(1i*(0:1:q));
dphi = my_ifft(testm,2*q+1);
% testm = phim*diag(1i*(0:1:q-1));
% dphi = my_ifft(testm,2*q);
%
% dphi/dr
drphim = 2*D*phim;
drphi = my_ifft(drphim,2*q+1);
% drphi = my_ifft(drphim,2*q);
%
% dq/dtheta
testm = qm*diag(1i*(0:1:q));
dq = my_ifft(testm,2*q+1);
% testm = qm*diag(1i*(0:1:q-1));
% dq = my_ifft(testm,2*q);
% dq/dr
drqm = 2*D*qm;
drq = my_ifft(drqm,2*q+1);
% drq = my_ifft(drqm,2*q);

% compute J_{q phi} = dQ/dr dphi/dtheta - dPhi/dr dq/dthete + dq/dr dphi/dtheta - dphi/dr dq/dtheta
Jqphi = diag(Q0_r)*dphi - diag(Phi0_r)*dq + Nonlinear*(drq.*dphi- drphi.*dq);


% dw/dtheta
testm = wm*diag(1i*(0:1:q));
dw = my_ifft(testm,2*q+1);
% testm = wm*diag(1i*(0:1:q-1));
% dw = my_ifft(testm,2*q);
%
% dw/dr
drwm = 2*D*wm;
drw = my_ifft(drwm,2*q+1);
% drw = my_ifft(drwm,2*q);

% compute J_{w phi} = dW/dr dphi/dtheta - dPhi/dr dw/dthete + dw/dr dphi/dtheta - dphi/dr dw/dtheta
Jwphi = diag(W0_r)*dphi - diag(Phi0_r)*dw + Nonlinear*(drw.*dphi-drphi.*dw);

% dpsi2/dtheta
testm = psi2m*diag(1i*(0:1:q));
dpsi2 = my_ifft(testm,2*q+1);
% dpsi2/dr
drpsi2m = 2*D*psi2m;
drpsi2 = my_ifft(drpsi2m,2*q+1);
%
% compute J_{psi2 q} = dPsi2/dr dq/dthete - dQ/dr dpsi2/dtheta + dpsi2/dr dq/dtheta -dq/dr dpsi2/dtheta  
Jpsi2q = diag(Psi20_r)*dq - diag(Q0_r)*dpsi2 + Nonlinear*(drpsi2.*dq-drq.*dpsi2);

% the Jacobian has a factor of 1/r in it; let's put it in.  We want to
% divide each row of the Jacobian matrix by 1/r.  This is done by
% premultiplying by a diagonal matrix whose diagonal is 1/riprime.
Jqphi = diag(1./riprime)*Jqphi;
Jwphi = diag(1./riprime)*Jwphi;
Jpsi2q = diag(1./riprime)*Jpsi2q;

% Now compute the FFT of the Jacobians.
Jqphi_m = my_fft(Jqphi,2*q+1);
Jwphi_m = my_fft(Jwphi,2*q+1);
Jpsi2q_m = my_fft(Jpsi2q,2*q+1);

% now throw away the K+1, ... 2K modes...
Jqphi_m = Jqphi_m(:,1:K+1);
Jwphi_m = Jwphi_m(:,1:K+1);
Jpsi2q_m = Jpsi2q_m(:,1:K+1);
% because I'm now truncating down from 0...K modes from 0...2K modes, I need
% to rescale the Fourier coefficients.
Jqphi_m = Jqphi_m*(2*K+1)/(2*q+1);
Jwphi_m = Jwphi_m*(2*K+1)/(2*q+1);
Jpsi2q_m = Jpsi2q_m*(2*K+1)/(2*q+1);


% for testing purposes...  build qold1 so that qold(1,:) = cos(theta).
% Then take qmold1 = my_fft(qold,2*K+1) and proceed here...
% testm = qm*diag(1i*(0:1:q));
% dq_phys = my_ifft(testm,2*q+1);
% the above should be -sin(theta)
% dq_phys(1,1)
% Jqphi = dq_phys + dq_phys.^2;
% Jqphi_m = my_fft(Jqphi,2*q+1);
% Jqphi_m = Jqphi_m(:,1:K+1);
% Jqphi_m = Jqphi_m*(2*K+1)/(2*q+1);
