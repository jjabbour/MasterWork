global tspan Nc K D2 D riprime dt Re Pr rr
init
% Initialize parameter needed for the time stepper 
% of the sheared electroconvection 
tspan = [0 .2];         % integration time  
K  = 32;                % highest Fourier wave
Nc  = 24;               % highest power of the Tchebyshev polynomial
rr = .56;               % aspect ratio
Re = .249;              % Reynolds number
Pr = 75.8;              % Prandlt number 
dt = 5.0e-4;            % time step 

% loading some precomputed input for the time stepper
% if the parameters (aspect ratio, Fourrier modes and Chebyschev) are
% changed we need to recompute the following output 
load Tq_vs_psi2_Nc24_K32.mat
load('base_state_charge_rrp56Nc24.mat')
InitBCs_amp  =  1.0000e-03;
[D,D2,yi,I,It,theta,riprime] = make_grids(Nc,K,rr);

% % if you are starting from a known solution leave these line commented
% % and uncomment the "load" line 
% % if you need to find an inital guess using long time integration 
% Ra = 580;  % initialize The initial dimensionless Rayleigh number
% [psi2_now,psi2m_now,q_now,qm_now,w_now,wm_now,phi_now,phim_now] = ...
%     random_initial_data(Tq_vs_psi2,D,D2,riprime,theta,Nc,K,InitBCs_amp);
  % or 
% % load('initial_guess_tspan_40_Ra_580_Re_p249');

% Reshaping the initial guess to meet the tailored time stepper and Cont_meth 
% U1 = reshape(psi2_now,(2*K+1)*(Nc+1),1); U1m = reshape(psi2m_now,(K+1)*(Nc+1),1);
% U2 = reshape(q_now,(2*K+1)*(Nc+1),1);    U2m = reshape(qm_now,(K+1)*(Nc+1),1);
% U3 = reshape(w_now,(2*K+1)*(Nc+1),1);    U3m = reshape(wm_now,(K+1)*(Nc+1),1);
% U4 = reshape(phi_now,(2*K+1)*(Nc+1),1);  U4m = reshape(phim_now,(K+1)*(Nc+1),1);
% U_init  = [U1;U2;U3;U4];
% Um_init = [U1m;U2m;U3m;U4m];

% time integration of the initial condition 
% [Uout,U_full] = Electro_time_stepper([0 10],U_init,Ra,dt,Pr,Re,rr,K,Nc);
% omega = .42;%(pi/6)/diff(tspan); % initializing a phase speed guess  

% plotting the output after integration
% % plotting the four physical quantities u = [q;\psi_2;\omega;\phi] after
% % integration 
% [psi2_f,q_f,w_f,phi_f] = vec_2_mat(Uout,Nc,K);
% [Psi20_phy,Q_phy,W_phy,Phi_phy] = vec_2_mat(U_full,Nc,K);

% figure()
% subplot(2,2,1) 
% polar3d(q_f,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
% title('charge q')
% subplot(2,2,2) 
% polar3d(psi2_f,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
% title('potential psi2')
% subplot(2,2,3) 
% polar3d(w_f,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
% title('vorticity w')
% subplot(2,2,4)
% polar3d(phi_f,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
% title('stream function phi') 
% figure()
% polar3d(Phi_phy,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
% pause

% starting from precomputed solution
load('tracing_periodic_solution_date735924.mat')
U_init = X(1:end-2,1);
omega  = X(end-1,1);
Ra     = X(end,1);

% Passign the arguments needed for the continuation solver
% imax is the number of solution point on the branch 
% ds is the step size in the parameter space 
% filename is the name of the file that stores the needded output 
imax = 9;
ds   = -1;
% The output file name; this file will be stored in the ..\mat_file 
filename = ['tracing_periodic_solution_ds_m1_try',int2str(today)];

[X,period,mu,V,eig_flag,bif_point] = CM_predictor_rel_equi([U_init;omega],Ra,filename,imax,ds);



