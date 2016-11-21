function [Psi2_norm,Q_norm,W_norm,Phi_norm] = norms(Psi2_phy,Q_phy,W_phy,Phi_phy,r,K,Nc);

dr =  - diff(r);
% these are the weights for the Trapezoid rule...  sum( delta_r.*f ) will
% approximate \int f(r) dr.
delta_r = [dr(1)/2; (dr(1:Nc-1)+dr(2:Nc))/2 ; dr(Nc)/2];
% sum(r.*delta_r)  % this approximates int r dr = r_out^2/2 - r_in^2/2 
A = pi*(r(1)^2 - r(Nc+1)^2);

dtheta = 2*pi/(2*K+1);

Psi2_norm = sqrt(sum(sum(dtheta*diag(r.*delta_r)*Psi2_phy.^2)));
Q_norm = sqrt(sum(sum(dtheta*diag(r.*delta_r)*Q_phy.^2)));
W_norm = sqrt(sum(sum(dtheta*diag(r.*delta_r)*W_phy.^2)));
Phi_norm = sqrt(sum(sum(dtheta*diag(r.*delta_r)*Phi_phy.^2)));



        