function [X,mu,V,eig_conv,bifpoint] = CM_predictor_fixed_point(X_init,Ra,filename,kmax,ds)
global tspan Pr Re rr Nc K 
% Continuation method  using natural parametrization 
% G(u,nu) = \phi(t,u(nu),nu) - u(nu) = 0, G: R^n x R -> R^n
% The nonlinear system G using Newton method and the linear system
% is solved using Arnoldi method through gmres
% D_u(G) * Delta_u = -G

% Input:
% u0 Nxm matrix with column vector initial condtion of the time stepper
% nu_init mx1 vector parameter value for respective initial condtion 
% Output:
% X N+1xm where X(:,m) =[u ; nu]; fixed point
if nargin < 3,    kmax = 50;  end
if nargin < 4,    ds = 1;     end
Nmax = 50;   % maximum Newton iteration
N = length(X_init); 
X = zeros(N+1,kmax);
mu = zeros(10,kmax);
V = zeros(N,10,kmax);
eig_conv = zeros(1,kmax);
time_count_eig = zeros(1,kmax);
NEW_res    = zeros(Nmax,kmax);
NEW_rel    = zeros(Nmax,kmax);
NEW_abs    = zeros(Nmax,kmax);
GMRES_rel  = zeros(Nmax,kmax);
GMRES_flag = zeros(Nmax,kmax);

bifpoint  = [];
run_info = ['The output for this file are based on the following parameters \n' ...
             'NEWTON tolerance:  abs_tol, res_tol and rel_tol: 1e-6 \n'...
             'GMRES residual tolerance 1e-3 and epsilon of the function evaluation is 1e-6\n'];

% initial guess for the rotation wave speed for saving purposes only          
% omega = X_init(end);

for k = 1:kmax
    [root,lambda,vec,flag_eigs,eig_time,New_res_vec,New_abs_vec,New_rel_vec,Gmres_rel_vec,Gmres_flag_vec]...
                     = New_Kry_solver_fix_point(X_init,Ra,Nmax);
    X(:,k)   = [root;Ra];               % solution vector [u,omega,Ra]
    mu(:,k)  = lambda;                  % eigenvalues of the map 
    V(:,:,k) = vec;                     % eigenvectors
    eig_conv(k) = flag_eigs;            % eigs flag checking for convergence of eigenvalues 
    time_count_eig(k)   = eig_time;     % CPU time for each computation 
    NEW_res(:,k)      = New_res_vec;    % Newton residual 
    NEW_abs(:,k)      = New_abs_vec;    % Newton absolute error
    NEW_rel(:,k)      = New_rel_vec;    % Newtom relative error
    GMRES_rel(:,k)    = Gmres_rel_vec;  % GMRES relative error
    GMRES_flag(:,k)   = Gmres_flag_vec; % GMRES flag
    
    % Initializng the guess for the new point on the solution branch
    Ra = Ra + ds;
    X_init = root;
    % checking for bifurcation point from stable known solution 
    if max(abs(lambda))> 1
        bifpoint = [bifpoint, k];
    end
    % saving data in a file 
    save(['..\mat_file\' filename],'X','mu','V','eig_conv','bifpoint','time_count_eig','NEW*','GMRES*',...
        'run_info','ds','kmax','tspan','Pr','Re','rr','Nc','K')
end
end

