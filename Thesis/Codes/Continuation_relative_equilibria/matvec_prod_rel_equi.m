function [Ax] = matvec_prod_rel_equi(unew,phi,u_theta_0,omega,Ra,Dv)
global Nc K tspan dt rr Re Pr 
tau = diff(tspan);
alpha = 0;
epsilon = 1e-6;
Dv1 = Dv(1:end-1);
Dv2 = Dv(end);

% forward finite difference for the Jacobian Matrix
vin  = unew +epsilon*Dv1;
phi_e_total = Electro_time_stepper(tspan,vin,Ra,dt,Pr,Re,rr,K,Nc);
phi_e = phi_e_total(:,end);
J = (phi_e-phi)/epsilon;

u_theta = dtheta_electro(unew,K,Nc);
gamma_deltau = rotation_symmetry(Dv1,omega,tspan,Nc,K);

Ax1 = J - gamma_deltau + tau*(u_theta*Dv2);
Ax2 = u_theta_0.'*Dv1 + alpha*Dv2;
Ax =[Ax1;Ax2];
end