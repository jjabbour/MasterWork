function [u_theta] = dtheta_electro_fft(u,K,Nc)
[psi2,q,w,phi] = vec_2_mat(u,Nc,K);

kwave = [0:K -K:-1];
A = zeros(Nc+1,length(kwave));
for k = 1:Nc+1
    A(k,:) = 1i*kwave;
end
A = A.';


psi2m = fft(psi2.');     psi2m_p = A.*psi2m;     psi2_p = ifft(psi2m_p).';
qm    = fft(q.');        qm_p    = A.*qm;        q_p    = ifft(qm_p).';
wm    = fft(w.');        wm_p    = A.*wm;        w_p    = ifft(wm_p).';
phim  = fft(phi.');      phim_p  = A.*phim;      phi_p  = ifft(phim_p).';

u_theta = mat_2_vec(psi2_p,q_p,w_p,phi_p,Nc,K);

end




