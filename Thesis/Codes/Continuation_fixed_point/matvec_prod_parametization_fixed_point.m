function [Ax] = matvec_prod_parametization_fixed_point(unew,phi,Ra,Dv)
global Nc K tspan dt rr Re Pr 
%tau = diff(tspan);
epsilon = 1e-6;

% forward finite difference for the Jacobian Matrix
vin  = unew +epsilon*Dv;
phi_e_total = Electro_time_stepper(tspan,vin,Ra,dt,Pr,Re,rr,K,Nc);
phi_e = phi_e_total(:,end);
J = (phi_e-phi)/epsilon;
Ax = J - Dv;

end