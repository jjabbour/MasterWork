function [psi2old2,psi2mold2,qold2,qmold2,phiold2,phimold2,wold2,wmold2,...
    Jqphi_hat2,Jwphi_hat2,Jpsi2q_hat2]=startup(psi2old1,psi2mold1,wold1,wmold1,...
    phiold1,phimold1,qold1,qmold1,Phi0_r,Q0_r,W0_r,Psi20_r,...
    riprime,Nonlinear,D,Nc,q);
% q = 2*K
psi2old2=psi2old1;
psi2mold2=psi2mold1;

qold2=qold1;
qmold2=qmold1;

phiold2=phiold1;
phimold2=phimold1;

wold2=wold1;
wmold2=wmold1;

% Compute the Jacobians at the past time-step...
[Jqphi_hat2,Jwphi_hat2,Jpsi2q_hat2] = compute_nonlinearity(psi2mold2,wmold2,...
    phimold2,qmold2,Phi0_r,Q0_r,W0_r,Psi20_r,riprime,Nonlinear,D,Nc,q);
