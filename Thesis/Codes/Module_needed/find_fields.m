function [Psi2_phy,psi2_phy,Psi2_r_phy,Q_phy,q_phy,Phi_phy,phi_phy,W_phy,w_phy,Ur_phy,Utheta_phy] ...
    = find_fields(Psi20_phy,Psi20_r_phy,Q0_phy,Phi0_phy,Phi0_r_phy,W0_phy,psi2m,qm,phim,wm,riprime,D,K)
% if pert = 1 then I get the full fields (including the base state).  Otherwise, 
% I just get the deviations from the base state
    
% r-component of the velocity is 1/r dphi/dtheta.
testm = phim*diag(1i*(0:1:K));
dphi = my_ifft(testm,2*K+1);
ur = diag(1./riprime)*dphi;
% theta-component of the velocity is -dphi/dr.  Factor of 2 is from the
% chain rule.
drphim = 2*D*phim;
drphi = my_ifft(drphim,2*K+1);
utheta = - drphi;
% recall that we've been computing deviations from the base state... so
% now we define the physical velocity
Ur_phy = ur; % Ur = 1/r dPhi0/dtheta + 1/r dphi/dtheta = 1/r dphi/dtheta
Utheta_phy = - Phi0_r_phy + utheta; % Utheta = - dPhi0/dr - dphi/dr

% construct the other quanties of interest
% the deviation of the potential from base state in physical space
psi2_phy = my_ifft(psi2m,2*K+1);  
% the deviation of the charge from base state in physical space
q_phy = my_ifft(qm,2*K+1); 
% the deviation of the stream function from base state in physical space
phi_phy = my_ifft(phim,2*K+1); 
% the deviation of the vorticity from base state in physical space
w_phy = my_ifft(wm,2*K+1);

% Now add these deviations to the base state to get the potential, charge, stream
% function, and vorticity in physical space.
Psi2_phy = psi2_phy + Psi20_phy;
Q_phy = q_phy + Q0_phy;
Phi_phy = phi_phy + Phi0_phy;
W_phy = w_phy + W0_phy;

% some things needed in order to compute the Nusselt number:
drpsi2m = 2*D*psi2m;
drpsi2 = my_ifft(drpsi2m,2*K+1);
Psi2_r_phy = drpsi2 + Psi20_r_phy;
end


