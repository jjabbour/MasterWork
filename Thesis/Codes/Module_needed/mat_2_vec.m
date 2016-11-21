function [Uout] =  mat_2_vec(psi2_phy,q_phy,w_phy,phi_phy,Nc,K)
 
U1 = reshape(psi2_phy,(Nc+1)*(2*K+1),1);
U2 = reshape(q_phy,(Nc+1)*(2*K+1),1);
U3 = reshape(w_phy,(Nc+1)*(2*K+1),1);
U4 = reshape(phi_phy,(Nc+1)*(2*K+1),1);
Uout = [U1;U2;U3;U4];
end
