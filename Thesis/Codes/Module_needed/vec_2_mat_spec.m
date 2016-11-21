function [psi2m_now,qm_now,wm_now,phim_now] =  vec_2_mat_spec(U_init,Nc,K)
len = (Nc+1)*(K+1);

psi2m_now = reshape(U_init(1:len),(Nc+1),(K+1));
qm_now    = reshape(U_init(len+1:2*len),(Nc+1),(K+1));
wm_now    = reshape(U_init(2*len+1:3*len),(Nc+1),(K+1));
phim_now  = reshape(U_init(3*len+1:4*len),(Nc+1),(K+1));

end