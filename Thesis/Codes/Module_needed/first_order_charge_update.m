function [psi2mnew,qmnew] = first_order_charge_update(psi2mnew,qmnew,qm_now,dt,Tq_vs_psi2,Lop,m,Nc,Jqphi_hat_now)


% at each wave number, want to take a time-step which is explicit in the
% nonlinearity and implicit in the diffusion term.
%
% (qm_new - qm_now)/dt + Jqphi_now = Lop psi2_new
%
% qm_new = qm_now - dt*Jqphi_now + dt*Lop*psi2_new
%
% qm_new - dt*Lop*psi2_new = qm_now - dt*Jqphi_now
%
% T psi2_new - dt*Lop psi2_new = qm_now - dt*Jqphi_now
%
% (T - dt*Lop) psi2_new = qm_now - dt*Jqphi_now.
%
%          -->  psi2new = inv( T - dt Lop) RHS
% We need to solve
%
%  A_full(2:Nc,1:Nc+1) psi2mnew(1:Nc+1,m+1) = B(2:Nc)
%
% subject to boundary conditions on psi2mnew(1,m+1) and psi2mnew(Nc+1,m+1)
%
% We want to solve Nc-1 equations (at the interior nodes)
% involving Nc+1 unknowns (including the boundary nodes).  
% The coefficients are given by the 2, ... Nc rows of A_full and the
% uknowns are psi2mnew from 1 to Nc+1.  But we know psi2mnew at 1 and Nc+1
% and so we transfer this to the RHS of the problem as follows:
psi2mnew(1,m+1)=0;
psi2mnew(Nc+1,m+1)=0;
A_full =  Tq_vs_psi2(:,:, m+1) - dt*Lop(:,:);
A_trunc = A_full(2:Nc,2:Nc);
B =  qm_now(2:Nc, m+1) - dt*Jqphi_hat_now(2:Nc, m+1);
RHS=B-[A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[psi2mnew(1,m+1);psi2mnew(Nc+1,m+1)];
psi2mnew(2:Nc,m+1) = A_trunc\RHS;
%
% Now, update the surface charge.
qmnew(:, m+1)=Tq_vs_psi2(:,:, m+1)*psi2mnew(:, m+1);


