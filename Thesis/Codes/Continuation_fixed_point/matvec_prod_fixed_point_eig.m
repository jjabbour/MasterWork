function [Av] = matvec_prod_fixed_point_eig(unew,phi,Ra,v)
% approximate the action of the Jacobian D_u\Phi_t(unew,Ra)*v 
% using a forward finite difference  

global tspan dt Pr Re rr K Nc
epsilon = 1e-6;
vin    = unew + epsilon*v;
phi_ue_total = Electro_time_stepper(tspan,vin,Ra,dt,Pr,Re,rr,K,Nc);
phi_ue = phi_ue_total(:,end);
Av = (phi_ue-phi)/epsilon;

end