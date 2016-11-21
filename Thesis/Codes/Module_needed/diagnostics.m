function [E_kin_dens,E_enstr_dens,Nu] ...
    = diagnostics(Ur_phy,Utheta_phy,W_phy,Q_phy,Psi20_r_phy,Psi2_r_phy,sigma,r,K,Nc)

dr =  - diff(r);
% these are the weights for the Trapezoid rule...  sum( delta_r.*f ) will
% approximate \int f(r) dr.
delta_r = [dr(1)/2; (dr(1:Nc-1)+dr(2:Nc))/2 ; dr(Nc)/2];
% sum(r.*delta_r)  % this approximates int r dr = r_out^2/2 - r_in^2/2 
A = pi*(r(1)^2 - r(Nc+1)^2);

dtheta = 2*pi/(2*K+1);
% theta = dtheta*(0:1:2*K);
% sum(cos(theta)*dtheta) % this approximates the integral of cos(theta)

E_kin = 1/2*(Ur_phy.^2 + Utheta_phy.^2); % Kinetic energy
E_enstr = 1/2*W_phy.^2;  % Enstrophy = 1/2* vorticity^2.

% E_kin = r*sin(theta).^2; for testing purposes; this converges to an integral of
% (7pi/3)/A with rate dtheta^2+dr^2
E_kin_dens = sum(sum(dtheta*diag(r.*delta_r)*E_kin))/A;
E_enstr_dens = sum(sum(dtheta*diag(r.*delta_r)*E_enstr))/A;

% sigma3 = 1.3*10^(-7); % from page 3624 of Phys Fluids 11(1999)
% one_layer_s = 3.6*10^(-9); % each layer is 3.6nm
% num_layers = 75; % from Stephen
% s = num_layers*one_layer_s;
% sigma = sigma3*s;
arg1 = -sigma*Psi2_r_phy + Q_phy.*Ur_phy;
arg2 = -sigma*Psi20_r_phy;
Nu = sum(dtheta*transpose(arg1))./sum(dtheta*transpose(arg2));


        