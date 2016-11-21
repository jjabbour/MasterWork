function [u_theta] = dtheta_electro(u,K,Nc)

kwave = 0:K;
A = zeros(Nc+1,length(kwave));
for k = 1:Nc+1
    A(k,:) = 1i*kwave;
end

[psi2,q,w,phi] = vec_2_mat(u,Nc,K);

psi2m = my_fft(psi2,2*K+1);     psi2m_p = A.*psi2m;     psi2_p = my_ifft(psi2m_p,2*K+1);
qm    = my_fft(q,2*K+1);        qm_p    = A.*qm;        q_p    = my_ifft(qm_p,2*K+1);
wm    = my_fft(w,2*K+1);        wm_p    = A.*wm;        w_p    = my_ifft(wm_p,2*K+1);
phim  = my_fft(phi,2*K+1);      phim_p  = A.*phim;      phi_p  = my_ifft(phim_p,2*K+1);

u_theta = mat_2_vec(psi2_p,q_p,w_p,phi_p,Nc,K);

end




