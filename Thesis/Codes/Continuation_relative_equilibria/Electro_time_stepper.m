function [U_perturbation,U_full] = Electro_time_stepper(tspan,U0,Ra,dt,Pr,Re,rr,K,Nc)
load Tq_vs_psi2_Nc24_K32.mat
load('base_state_charge_rrp56Nc24.mat')

if nargin < 4,    dt  = 5.0e-04;    end  % size of the time step:
if nargin < 5,    Pr  = 75.8;       end  % Prandtl number: 
if nargin < 6,    Re  = 0.249;      end  % Reynolds number:
if nargin < 7,    rr  = .56;        end  % relative ratio
if nargin < 8,    K   = 32;         end  % number of Fourier modes in the angular direction
if nargin < 9,    Nc  = 24;         end  % number of meshpoints in the radial direction:


Omega = Re*Pr*(1-rr)/rr;    % the rotation rate of the inner electrode.  Need this to compute the base state.
Nrun = round((tspan(2)-tspan(1))/dt);
% The following are used by the script diagnostics.m which computes diagnostic quantities such 
% as the kinetic energy density, the enstrophy density, and the Nusselt number all as a function of time.
% charge relaxation timescale
eps0        = 8.85*10^(-12);    %permittivity of free space
sigma3      = 1.3*10^(-7);      % from page 3624 of Daya Phys Fluids 11(1999)
one_layer_s = 3.6*10^(-9);      % each layer is 3.6nm
num_layers  = 75;               % from Stephen
s           = num_layers*one_layer_s;
sigma       = sigma3*s;
r_out       = 6.4*10^(-3);      % from page 3624 of Daya Phys Fluids 11(1999)
r_in        = rr*r_out;         % by definition of the radius ratio rr
d           = r_out - r_in;     % by definition of d
tau_q       = eps0*d/sigma;     % the charge relaxation time

% In the time-stepping, we need to update the vorticity.  This uses BCs4dr_phi
BCs4dr_phi  = [0;0];
% to turn off the nonlinear terms for diagnostic reasons, set Nonlinear to zero.
Nonlinear   = 1;

% ==============================================================
% Create Grids and the differentiation Matrices
% ==============================================================
[D,D2,yi,I,It,theta,riprime] = make_grids(Nc,K,rr);
% ==============================================================
% Instead of loading passed as global variable
% load Tq_vs_psi2_Nc24_K32.mat
% compute the base state charge and potential using the symbolic toolbox --- this is slow!
% [Psi20,Psi20_r,Q0,Q0_r] = base_state_charge(rr, riprime); 
% or load them directly from a file
% load base_state_charge_Nc24.mat
%load('base_state_charge_rrp56Nc24.mat')
[Phi0,Phi0_r,Phi0_rr,W0,W0_r] = base_state_flow(Omega, rr, riprime); 

% As given, the base state quantities are independent of theta; they're
% (Nc+1) x 1 matrices.  We want to create (Nc+1) x (2K+1) matrices which
% will represent the base state quantities in the annulus for use in
% reconstructing the unknowns in the annulus (find_fields.m) as well as in
% computing diagnostic quantities (diagnostics.m and norms.m)
Psi20_phy = Psi20*ones(1,2*K+1);
Psi20_r_phy = Psi20_r*ones(1,2*K+1);
Q0_phy = Q0*ones(1,2*K+1);
Phi0_phy = Phi0*ones(1,2*K+1);
Phi0_r_phy = Phi0_r*ones(1,2*K+1);
W0_phy = W0*ones(1,2*K+1);

% =======================================================
%   Initial Conditions and start-up
% =======================================================
% initialize arrays that would otherwise dynamically grow:
% [Ur_phy,Utheta_phy,phi_phy,psi2_phy,q_phy,w_phy,time,E_kin,E_enstr,Nu,psi2_norm,q_norm,phi_norm,w_norm] = initialize(Nc,K,Nrun);
% for random initial data:
% random_initial_data_PT
%[psi2_now,psi2m_now,q_now,qm_now,w_now,wm_now,phi_now,phim_now] = random_initial_data(Tq_vs_psi2,D,D2,riprime,theta,Nc,K,InitBCs_amp);
[psi2m_now,qm_now,wm_now,phim_now] =  reshape_IC_time_stepper(U0,Nc,K);
%[psi2_now,q_now,w_now,phi_now] =  vec_2_mat(U0,Nc,K);
% if you want to use the same initial data for multiple runs then load it
% from a file:
% load rand_ID.mat
% t_start = 0;
% If you want to continue a run from data that you saved already in a file
% [Ra,Pr,Re,rr,psi2_now,psi2m_now,q_now,qm_now,phi_now,phim_now,w_now,wm_now,t_start] = load_ID('Pr10_Ra70_Re0p2_rr0p56.mat');

% Why am I using my own fft and ifft you may well ask?  Well the arrays in
% this code are being stored in the opposite way that matlab's fft/ifft are
% expecting.  And so I have to do transposes and stuff like that.  This is
% a vestige of how Peichun wrote the code; I think she stored the data in
% this particular way so that she could use polar3d.m to plot the data.

% in case you want to see the spectra of the initial data
% plot_spectra
% pause(1)

% precompute the influence matrix, needed for the vorticity/fluid potential
% solve.
[wm1,phim1,wm3,phim3,Mi] = influence_matrix(D,D2,riprime,K,dt,Pr,It);

fprintf('//Ra = %.2f:\n', Ra);
fprintf('//Pr = %.2f:\n', Pr);
fprintf('//Re = %.2f:\n', Re);

% ===================== Begin the time iteration ============================
% we're using a three-level timestepping scheme.  So we start with one
% step of a first-order accurate two-level scheme.  And then continue on
% with the three-level scheme.

% take one step using first-order time-stepping
% need to initialize these variables.
psi2m_new = zeros(Nc+1,K+1);
qm_new = zeros(Nc+1,K+1);
phim_new = zeros(Nc+1,K+1);
wm_new = zeros(Nc+1,K+1);
% compute the initial terms that go into the nonlinearity.
[Jqphi_hat_now,Jwphi_hat_now,Jpsi2q_hat_now] = compute_nonlinearity(psi2m_now,wm_now,...
    phim_now,qm_now,Phi0_r,Q0_r,W0_r,Psi20_r,riprime,Nonlinear,D,Nc,2*K);
[psi2m_new,qm_new,phim_new,wm_new] = first_order_time_step(psi2m_new,...
    qm_now,qm_new,phim_new,wm_now,wm_new,Jqphi_hat_now,...
    Jwphi_hat_now,Jpsi2q_hat_now,D,D2,K,riprime,...
    (dt)^2,Tq_vs_psi2,Nc,BCs4dr_phi,wm1,wm3,phim1,phim3,Mi,Ra,Pr,It);
% update so we can take the next timestep.
psi2m_old = psi2m_now;
qm_old = qm_now;
wm_old = wm_now;
phim_old = phim_now;
Jqphi_hat_old = Jqphi_hat_now;
Jwphi_hat_old = Jwphi_hat_now;
Jpsi2q_hat_old = Jpsi2q_hat_now;
% to do the not-smart option of no first-order time-step, comment out the
% following lines...
psi2m_now = psi2m_new;
qm_now = qm_new;
wm_now = wm_new;
phim_now = phim_new;

for jj=2:Nrun
    % compute the Jacobian terms pseudospectrally
    [Jqphi_hat_now,Jwphi_hat_now,Jpsi2q_hat_now] = compute_nonlinearity(psi2m_now,wm_now,...
        phim_now,qm_now,Phi0_r,Q0_r,W0_r,Psi20_r,riprime,Nonlinear,D,Nc,2*K);
    % Now timestep in Fourier space...
    for m=0:1:K % Loop for each Fourier mode
        % Lop_PT = 4*D2 + 2*riprime_invM.*D - m^2* full(sparse(1:Nc+1, 1:Nc+1, 1./riprime_square )); % Laplace operator matrix (checked)
        Lop = 4*D2 + 2*diag(1./riprime)*D - m^2*diag(1./(riprime.^2));
        % first step: solve for the new surface charge and the new 2d
        % electric potential
        [psi2m_new,qm_new] = second_order_charge_update(psi2m_new,qm_new,qm_now,qm_old,dt,...
            Tq_vs_psi2,Lop,m,Nc,Jqphi_hat_now,Jqphi_hat_old);
        % second step: now solve for the new vorticity and fluid potential
        [wm_new,phim_new] = second_order_vorticity_update(wm_new,phim_new,wm_now,wm_old,D,...
            BCs4dr_phi,wm1,wm3,phim1,phim3,Mi,Ra,Pr,dt,It,Lop,m,Nc,...
            Jwphi_hat_now,Jwphi_hat_old,Jpsi2q_hat_now,Jpsi2q_hat_old);
    end
    % update so we can take the next timestep.
    psi2m_old = psi2m_now;
    qm_old = qm_now;
    wm_old = wm_now;
    phim_old = phim_now;
    Jqphi_hat_old = Jqphi_hat_now;
    Jwphi_hat_old = Jwphi_hat_now;
    Jpsi2q_hat_old = Jpsi2q_hat_now;
    psi2m_now = psi2m_new;
    qm_now = qm_new;
    wm_now = wm_new;
    phim_now = phim_new;
end
%ii = 1;
% save data profiles into the arrays that we'll ultimately save in a
% data file...
[Psi2_phy,psi2_phy,Psi2_r_phy,Q_phy,q_phy,Phi_phy,phi_phy,W_phy,w_phy,Ur_phy,Utheta_phy] ...
        = find_fields(Psi20_phy,Psi20_r_phy,Q0_phy,Phi0_phy,Phi0_r_phy,W0_phy,psi2m_now,qm_now,phim_now,wm_now,riprime,D,K);
% [Psi2_phy,psi2_phy(:,:,ii+1),Psi2_r_phy,Q_phy,q_phy(:,:,ii+1),Phi_phy(:,:,ii+1),phi_phy(:,:,ii+1),W_phy,w_phy(:,:,ii+1),Ur_phy(:,:,ii+1),Utheta_phy(:,:,ii+1)] ...
%                               = find_fields(Psi20_phy,Psi20_r_phy,Q0_phy,Phi0_phy,Phi0_r_phy,W0_phy,psi2m_now,qm_now,phim_now,wm_now,riprime,D,K);
% [E_kin(ii+1),E_enstr(ii+1),Nu(:,ii+1)] = diagnostics(Ur_phy(:,:,ii+1),Utheta_phy(:,:,ii+1),W_phy,Q_phy,Psi20_r_phy,Psi2_r_phy,sigma,riprime,K,Nc);
% [psi2_norm(ii+1),q_norm(ii+1),w_norm(ii+1),phi_norm(ii+1)] = norms(psi2_phy(:,:,ii+1),q_phy(:,:,ii+1),W_phy,Phi_phy(:,:,ii+1),riprime,K,Nc);

U_full = mat_2_vec(Psi2_phy,Q_phy,W_phy,Phi_phy,Nc,K);
U_perturbation = mat_2_vec(psi2_phy,q_phy,w_phy,phi_phy,Nc,K);
end


function [psi2m_last,qm_last,wm_last,phim_last,Um_init] = reshape_IC_time_stepper(unew,Nc,K)
% [psi2,q,w,phi] 2D potential , charge, vorticity, stream line, each (2K+1)x(Nc+1)
% Input 
% unew  4(2K+1)(Nc+1) x  1   
% Nc    1 x 1                   Highest order of Chebychev polynomial
% K     1 x 1                   Highest wave number of Fourrier bases function 

% Output:
% Um_init   4(K+1)(Nc+1) x 1    Spectral transform of unew to be initial conditiion for the time stepper implemented 

m = length(unew);
q = 4*(2*K+1)*(Nc+1);
if m ~= q
    unew = unew(1:q,1);
end
    
[psi2_last,q_last,w_last,phi_last] = vec_2_mat(unew,Nc,K);
psi2m_last = my_fft(psi2_last,2*K+1);   U1m = reshape(psi2m_last,(K+1)*(Nc+1),1);
qm_last = my_fft(q_last,2*K+1);         U2m = reshape(qm_last,(K+1)*(Nc+1),1);
wm_last = my_fft(w_last,2*K+1);         U3m = reshape(wm_last,(K+1)*(Nc+1),1);
phim_last = my_fft(phi_last,2*K+1);     U4m = reshape(phim_last,(K+1)*(Nc+1),1);
Um_init = [U1m;U2m;U3m;U4m];

end