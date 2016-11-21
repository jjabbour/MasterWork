function [gamma] = rotation_symmetry(u,omega,tspan,Nc,K)

kwave = 0:K;
tau = diff(tspan);
A = zeros(Nc+1,length(kwave));
for k = 1:Nc+1
    A(k,:) = 1i*kwave*omega*tau;
end
B = exp(-A);

[psi2,q,w,phi] = vec_2_mat(u,Nc,K);
psi2m = my_fft(psi2,2*K+1);
qm    = my_fft(q,2*K+1);
wm    = my_fft(w,2*K+1);
phim  = my_fft(phi,2*K+1);

psi2m_rot = B.*psi2m;
qm_rot    = B.*qm;
wm_rot    = B.*wm;
phim_rot  = B.*phim;

psi2_now = my_ifft(psi2m_rot,2*K+1);
q_now    = my_ifft(qm_rot,2*K+1);
w_now    = my_ifft(wm_rot,2*K+1);
phi_now  = my_ifft(phim_rot,2*K+1);

gamma = mat_2_vec(psi2_now,q_now,w_now,phi_now,Nc,K);
end



