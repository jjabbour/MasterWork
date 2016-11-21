% This function returns the base state stream function and vorticity for a constant shear rate, Omega.
function [Phi0,Phi0_r,Phi0_rr,W0,W0_r] = base_state_flow(Omega, alpha, r)

% Inputs: Omega - dimensionless shear rate, alpha - radius ratio, r - the grid points
A = alpha^2*Omega/(1-alpha^2);

% Note that the base state isn't zero at the inner and outer radii.
Phi0 = A*(r.^2/2 - 1/(1-alpha)^2*log(r));
% set the potential so that Phi0 = 0 at the outer electrode.  This means
% that the total flux along a segment connecting the inner to the outer
% electrode will be Phi0 evaluated at the inner electrode...
Phi0 = Phi0 - Phi0(1); %*ones(size(Phi0));
Phi0_r = A*(r - 1/(1-alpha)^2*1./r);
Phi0_rr = A*(1 + 1/(1-alpha)^2*1./r.^2);

W0 = -2*A*ones(size(Phi0));
% Here is the radial derivative of the base state vorticity.  It's here
% just in case we ever want to perturb off of something other than the base
% state.
W0_r = zeros(size(Phi0));
end
