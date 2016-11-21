function [wmnew,phimnew] = second_order_vorticity_update(wmnew,phimnew,wm_now,wm_old,D,BCs4dr_phi,wm1,wm3,phim1,phim3,Mi,Ra,Pr,dt,It,Lop,m,Nc,Jwphi_hat_now,Jwphi_hat_old,Jpsi2q_hat_now,Jpsi2q_hat_old)

% first step: solve
%           (I-delta L) w = RHS
%                    L  phi = -w
% subject to the boundary conditions w(\pm 1)=0 and phi(\pm 1) = 0.
%
% Really, we'd like the boundary conditions to be either
%          w(\pm 1) = 0 and phi(1) = 0 and phi(-1) = g2
% or we'd like them to be
%          w(\pm 1) = 0 and phi(1) = 0 and phi_theta(-1) = 0
A_full = (It - 2/3*dt*Pr*Lop(:,:));
A_trunc = A_full(2:Nc,2:Nc);
% here are the boundary conditions on the vorticity
wm2(1,m+1)=0;
wm2(Nc+1,m+1)=0;
% We use these to adjust the RHS...
B = 4/3*wm_now(2:Nc, m+1) - 1/3*wm_old(2:Nc, m+1)...
    - 2/3*dt*(2*Jwphi_hat_now(2:Nc, m+1) - Jwphi_hat_old(2:Nc, m+1) )...
    + 2/3*dt*Pr*Ra*(2*Jpsi2q_hat_now(2:Nc, m+1) - Jpsi2q_hat_old(2:Nc, m+1)) ;
RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[wm2(1,m+1);wm2(Nc+1,m+1)];
wm2(2:Nc,m+1) = A_trunc\RHS;
% Now compute the potential from wm2...
A_full = Lop;
A_trunc = A_full(2:Nc,2:Nc);
% here are the boundary conditions on the potential
phim2(1,m+1)=0;
phim2(Nc+1,m+1)=0; % this is where we'd put the g2
B = -wm2(2:Nc,m+1);
RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[phim2(1,m+1);phim2(Nc+1,m+1)];
phim2(2:Nc,m+1) = A_trunc\RHS;

if norm(phim2) ==0
    % display('hello')
    BCs4w(:,:,m+1) = [0; 0];
else
    % I don't understand why BCs4dr_phi is here....
    BCs4w (:,:, m+1) = Mi(:,:,m+1)\( BCs4dr_phi - [2*D(Nc+1, :)*phim2(:,m+1); 2*D(1, :)*phim2(:, m+1)] );
    % BCs4w_PT(:,:, m+1) = Mi(:,:,m+1)\( BCs4dr_phi - [2*D(Nc+1, :)*phim2_PT(:,m+1); 2*D(1, :)*phim2_PT(:, m+1)] );
end

wmnew(:,m+1)=wm2(:,m+1)+BCs4w(1,1,m+1)*wm1(:,m+1)+BCs4w(2,1,m+1)*wm3(:,m+1);
phimnew(:,m+1)=phim2(:,m+1)+BCs4w(1,1,m+1)*phim1(:,m+1)+BCs4w(2,1,m+1)*phim3(:,m+1);

