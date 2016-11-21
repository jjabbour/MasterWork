function [psi2_now,q_now,w_now,phi_now] =  vec_2_mat(U_init,Nc,K)
len = (Nc+1)*(2*K+1);

psi2_now = reshape(U_init(1:len),(Nc+1),(2*K+1));
q_now    = reshape(U_init(len+1:2*len),(Nc+1),(2*K+1));
w_now    = reshape(U_init(2*len+1:3*len),(Nc+1),(2*K+1));
phi_now  = reshape(U_init(3*len+1:4*len),(Nc+1),(2*K+1));

end