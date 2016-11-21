load('random_IC_K32Nc_24.mat')
K = 32; Nc = 24;
n = (2*K+1)*(Nc+1); 
len_u = 4*n;
[psi2m_now,qm_now,wm_now,phim_now] =  reshape_IC_time_stepper(U_init,Nc,K);
Um_init = mat_2_vec_spec(psi2m_now,qm_now,wm_now,phim_now,Nc,K);
psi2_norm = norm(U_init(1:n));
q_norm = norm(U_init(n+1:2*n));
w_norm = norm(U_init(2*n+1:3*n));
phi_norm = norm(U_init(3*n+1:4*n));

nm = (K+1)*(Nc+1);
psi2m_norm = norm(Um_init(1:nm));
qm_norm = norm(Um_init(nm+1:2*nm));
wm_norm = norm(Um_init(2*nm+1:3*nm));
phim_norm = norm(Um_init(3*nm+1:4*nm));

NORM_Real = [psi2_norm;q_norm;w_norm;phi_norm];
NORM_Spec = [psi2m_norm;qm_norm;wm_norm;phim_norm];
test = NORM_Real - 2*NORM_Spec
