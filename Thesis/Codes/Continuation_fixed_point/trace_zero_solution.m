% clear; close all;clc;
init    % initializre all the modules needed in different folders
global tspan Nc K D2 D riprime dt Re Pr rr
% module needed for the time stepper
load Tq_vs_psi2_Nc24_K32.mat
load('base_state_charge_rrp56Nc24.mat')

% parameter needed for the time stepper of the sheared electroconvection 
tspan = [0 .2];         % time of integration 
K  = 32;                % highest fourier wave
Nc  = 24;               % highest power of the Chebyshev
rr = .56;               % aspect ratio
Re = .249;              % Dimensionless ratio number
Pr = 75.8;              % Dimensionless Prandlt number 
dt = 5.0e-4;            % time step 
InitBCs_amp  =  1.0000e-03;
               
[D,D2,yi,I,It,theta,riprime] = make_grids(Nc,K,rr);

Ra = 550;        % initialize The initial dimensionless Rayleigh number

% initializing random initial condition to time step to obtain 
% a good approximation for an initial solution for the continuation method  
% if you want to start from a known solution comment the lines that follow
[psi2_now,psi2m_now,q_now,qm_now,w_now,wm_now,phi_now,phim_now] = ...
    random_initial_data(Tq_vs_psi2,D,D2,riprime,theta,Nc,K,InitBCs_amp);
U1 = reshape(psi2_now,(2*K+1)*(Nc+1),1); U1m = reshape(psi2m_now,(K+1)*(Nc+1),1);
U2 = reshape(q_now,(2*K+1)*(Nc+1),1);    U2m = reshape(qm_now,(K+1)*(Nc+1),1);
U3 = reshape(w_now,(2*K+1)*(Nc+1),1);    U3m = reshape(wm_now,(K+1)*(Nc+1),1);
U4 = reshape(phi_now,(2*K+1)*(Nc+1),1);  U4m = reshape(phim_now,(K+1)*(Nc+1),1);
U_init  = [U1;U2;U3;U4];
Um_init = [U1m;U2m;U3m;U4m];
[Uout,U_full] = Electro_time_stepper([0 4],U_init,Ra,dt,Pr,Re,rr,K,Nc);
[psi2_f,q_f,w_f,phi_f] = vec_2_mat(Uout,Nc,K);
[Psi20_phy,Q_phy,W_phy,Phi_phy] = vec_2_mat(U_full,Nc,K);

% plotting the four physical quantities u = [q;\psi_2;\omega;\phi] after
% intgrration 
figure()
subplot(2,2,1) 
polar3d(q_f,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
title('charge q')
subplot(2,2,2) 
polar3d(psi2_f,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
title('potential psi2')
subplot(2,2,3) 
polar3d(w_f,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
title('vorticity w')
subplot(2,2,4)
polar3d(phi_f,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
title('stream function phi')

figure(2)
polar3d(Phi_phy,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
pause

% Passign the arguments needed for the continuation solver
% kmax is the number of solution point on the branch 
% ds is the step size in the parameter space 
% filename is the name of the file that stores the needded output if you
% want see the output received check CM_predictor
kmax = 130;
ds = 1;
filename = ['following_zero_solution_branch',int2str(today)];

[X,mu,V,eig_flag,bif_point] = CM_predictor_fixed_point(Uout,Ra,filename,kmax,ds);
% [X,mu,V,eig_flag,bif_point] = CM_predictor(X(1:end-1,26),X(end,26),filename,kmax,ds);