function [Av] = matvec_prod_rel_equi_eig(u_star,phi_period,Ra,period,v)
global dt Pr Re rr K Nc
epsilon = 1e-6;

vin    = u_star + epsilon*v;
phi_ue_total = Electro_time_stepper([0, period],vin,Ra,dt,Pr,Re,rr,K,Nc);
phi_ue = phi_ue_total(:,end);
Av = (phi_ue-phi_period)/epsilon;

end