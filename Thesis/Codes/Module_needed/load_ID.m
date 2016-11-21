function [Ra,Pr,Re,rr,psi2_now,psi2m_now,q_now,qm_now,phi_now,phim_now,w_now,wm_now,t_start] = load_ID(file_name)

load(file_name)

[a,b,c] = size(w_phy);
% a = Nc+1
% b = 2*K+1
% c = number of solutions saved; want the last one

psi2_now = psi2_phy(:,:,c);
q_now = q_phy(:,:,c);
phi_now = phi_phy(:,:,c);
w_now = w_phy(:,:,c);

psi2m_now = my_fft(psi2_now,b);
qm_now = my_fft(q_now,b);
phim_now = my_fft(phi_now,b);
wm_now = my_fft(w_now,b);

t_start = max(time);

