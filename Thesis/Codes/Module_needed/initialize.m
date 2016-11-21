function [Ur_phy,Utheta_phy,phi_phy,psi2_phy,q_phy,w_phy,time,E_kin,E_enstr,Nu,psi2_norm,q_norm,phi_norm,w_norm] = initialize(Nc,K,Nrun)

Ur_phy = zeros(Nc+1,2*K+1,Nrun+1);
Utheta_phy = zeros(Nc+1,2*K+1,Nrun+1);
phi_phy = zeros(Nc+1,2*K+1,Nrun+1);
psi2_phy = zeros(Nc+1,2*K+1,Nrun+1);
q_phy = zeros(Nc+1,2*K+1,Nrun+1);
w_phy = zeros(Nc+1,2*K+1,Nrun+1);

time = zeros(1,Nrun+1);
E_kin = zeros(1,Nrun+1);
E_enstr = zeros(1,Nrun+1);
Nu = zeros(Nc+1,Nrun+1);

psi2_norm = zeros(1,Nrun+1);
q_norm = zeros(1,Nrun+1);
phi_norm = zeros(1,Nrun+1);
w_norm = zeros(1,Nrun+1);
