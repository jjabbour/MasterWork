init
load('tracing_periodic_solution_date735924')
bif_point = 57; % bifurcation to amplitude vascillating occurs in column 57 of the results 
a = 0; % number of the Rayleigh number chosen from the output file (a+1) 
tf = 30;
run_range = bif_point+2; 
% run_range = bif_point:bif_point+a; %
Ra_v   = X(end,run_range);  % Collecting the Rayleigh number  corresponding to the run needed

% obtaining the unstable periodic solution computed by Cont_Met
U_star = X(1:end-2,65);     
% initialize the initial condition by adding a random initial condition to
% the unstable periodic solution
U_init = U_star + rand(size(U_star))*1e-2; 
% initialize the mat file name where the output will be stored
filename = ['direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_',int2str(tf),'Ra_',int2str(min(Ra_v)),'to',int2str(max(Ra_v))];
% filename = 'test';

% initializing the dimension of sim_output
m = size(U_star,1);
n = 5/5e-4;
p = length(Ra_v);
% 3D object 'sim_output' size is  (row,col,z) = (4*(Nc+1)*(2*K+1) , tf/dt, a+1) dt = 5e-4;
% corresponding to the solution U integrated over time tf for the
% respective number Ra from Ra_v
sim_output = zeros(m,n,p);
% running the simulation
for k = 1:p
    U_t = Electro_periodic_check([0,tf],U_init,Ra_v(k));
    sim_output(:,:,k) = U_t; % saving all the time steps 
    sim_output(:,:,k) = U_t(:,end-n+1:end); % saving the last 10000 time steps
    save(filename, 'U_*','sim_output','tf','Ra_v','Re','Pr','rr','K','Nc','-v7.3')
end

