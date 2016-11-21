% This function returns the base state stream function and vorticity for a constant shear rate, Omega.
function [Psi20,Psi20_r,Q0,Q0_r] = base_state_charge(alpha, r)

% Inputs: alpha - radius ratio, r - the grid points

Psi20 = (log(1-alpha)+log(r))/log(alpha);
% My Psi20_r agrees with the one Peichun had stored in Psi20_r.mat
Psi20_r = 1./(r*log(alpha));

r_out = r(1);
Nc = length(r)-1;
r_in = r(Nc+1);
a = 1/2;
b = 1/2;
c = 1;
% Preallocate for speed 
Q0 = zeros(Nc+1,1); Q0_r = zeros(Nc+1,1);
Q0(2:Nc,1) = (1/r_out).*hypergeom([a,b],c,r(2:Nc).^2/r_out^2) ...
          - (1./r(2:Nc)).*hypergeom([a,b],c,r_in^2./r(2:Nc).^2);
Q0 = Q0*(2/log(alpha));
Q0(Nc+1) = 0;

Q0_r(2:Nc,1) = (1/r_out)*(a*b/c).*hypergeom([a+1,b+1],c+1,r(2:Nc).^2/r_out^2).*(2*r(2:Nc)/r_out^2)...
    - (-1./r(2:Nc).^2).*hypergeom([a,b],c,r_in^2./r(2:Nc).^2)...
    - (1./r(2:Nc)).*((a*b/c)*hypergeom([a+1,b+1],c+1,r_in^2./r(2:Nc).^2)).*(-2*r_in^2./r(2:Nc).^3);
Q0_r(Nc+1)=0;
Q0_r = Q0_r*(2/log(alpha));

end
