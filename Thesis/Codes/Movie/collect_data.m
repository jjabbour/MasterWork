init
load('tracing_periodic_solution_date735924')
x_k = [49 53 61 65];
n = length(x_k);
X_star = X(:,x_k);
U_init = X_star(1:end-2,:);
Ra = X_star(end,:);
% the file name of the file that the output will be stored 
filename = ['time_stepping1em1_Ra_' int2str(Ra)];

tspan = [0 10];         % time of integration 
K  = 32;                % highest Fourrier wave
Nc  = 24;               % highest power of the Chebyshev
rr = .56;               % aspect ratio
Re = .249;              % Dimensionless ratio number
Pr = 75.8;              % Dimensionless Prandlt number 
dt = 5.0e-4;            % time step 
U_init_per = U_init + rand(size(U_init))*1e-1;
U_per   = zeros(4*(Nc+1)*(2*K+1),round(diff(tspan)/dt),n);

for k = 1:n
    sol = Electro_periodic_check(tspan,U_init_per(:,k),Ra(k),dt,Pr,Re,rr,K,Nc);
    U_per(:,:,k) = sol;
end
save(filename,'X_star','U*','rr','Ra','Re','dt','K','Nc','Pr','-v7.3')
