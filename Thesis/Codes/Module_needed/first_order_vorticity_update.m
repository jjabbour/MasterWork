function [wmnew,phimnew] = first_order_vorticity_update(wmnew,phimnew,wm_now,D,BCs4dr_phi,wm1,wm3,phim1,phim3,Mi,Ra,Pr,dt,It,Lop,m,Nc,Jwphi_hat_now,Jpsi2q_hat_now)

% One step scheme, doing the nonlinear term explicitly and the diffusion
% term implicitly
%
% (w_new-w_now)/dt + Jphiw_now = Pr L w_new + Pr Ra Jpsi2_now
%
% w_new-w_now + dt Jphiw_now = dt Pr L w_new + dtPr Ra Jpsi2_now
%
% w_new - dt Pr L w_new = w_now - dt Jphiw_now  + dt Pr Ra Jpsi2_now
%
% (I - dt Pr L ) w_new = w_now - dt Jphiw_now  + dt Pr Ra Jpsi2_now
%
% first step: solve
%           (I-delta L) w = RHS
%                    L  phi = -w
% subject to the boundary conditions w(\pm 1)=0 and phi(\pm 1) = 0.
%
% Really, we'd like the boundary conditions to be either
%          w(\pm 1) = 0 and phi(1) = 0 and phi(-1) = g2
% or we'd like them to be
%          w(\pm 1) = 0 and phi(1) = 0 and phi_theta(-1) = 0
A_full = (It - dt*Pr*Lop(:,:));
A_trunc = A_full(2:Nc,2:Nc);
% here are the boundary conditions on the vorticity
wm2(1,m+1)=0;
wm2(Nc+1,m+1)=0;
% We use these to adjust the RHS...
B = wm_now(2:Nc, m+1) - dt*Jwphi_hat_now(2:Nc, m+1)...
    + dt*Pr*Ra*Jpsi2q_hat_now(2:Nc, m+1);
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

