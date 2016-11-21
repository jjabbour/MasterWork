function [psi2mnew,qmnew] = second_order_charge_update(psi2mnew,qmnew,qm_now,qm_old,dt,Tq_vs_psi2,Lop,m,Nc,Jqphi_hat_now,Jqphi_hat_old)


% at each wave number, want to solve
%   q_new - 2dt/3 Lop psi2new = 4/3 q_now - 1/3 q_old 
%                                      - 4dt/3 Jqphi_now + 2dt/3 Jqphi_old
% Recalling that q = T psi2, this becomes
%   T psi2new - 2dt/3 Lop psi2new = RHS
%          -->  psi2new = inv( T - 2dt/3 Lop) RHS
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
A_full =  Tq_vs_psi2(:,:, m+1) - 2/3*dt*Lop(:,:);
A_trunc = A_full(2:Nc,2:Nc);
B =  4/3*qm_now(2:Nc, m+1) - 1/3*qm_old(2:Nc, m+1)...
            - 2/3*dt*(2*Jqphi_hat_now(2:Nc, m+1) - Jqphi_hat_old(2:Nc, m+1)) ;
RHS=B-[A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[psi2mnew(1,m+1);psi2mnew(Nc+1,m+1)];
psi2mnew(2:Nc,m+1) = A_trunc\RHS;
%
% Now, update the surface charge.
qmnew(:, m+1)=Tq_vs_psi2(:,:, m+1)*psi2mnew(:, m+1);

% To look at the condition numbers of the matrices we're inverting...      
% check_condition_numbers
