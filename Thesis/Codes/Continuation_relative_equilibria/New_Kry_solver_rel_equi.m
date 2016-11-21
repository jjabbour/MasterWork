function [root,period,mu,V,flag_eigs,eig_time,New_res_vec,New_abs_vec,New_rel_vec,Gmres_rel_vec,Gmres_flag_vec] ...
                        = New_Kry_solver_rel_equi(X0,Ra,imax,abs_tol,rel_tol)
% MVNewton Newton method using gmres to find the root iteratively of the nonlinear equation 
% [root,mu,New_iter,abs_err,rel_err] = New_Kry_solver_parametization(X0,imax,abs_tol,rel_tol)

% INPUT:
% X0    mx1 initial point 
% abs_tol   tolerance of the error set at 1e-6 by default
% imax  number of iteration set at 50 by default
% OUTPUT:
% root  1xm root of the function
% iter     1x1 number of itteration
% epsi  1x1 the relative error 
global tspan Nc K dt Re Pr rr 
format shorte
if nargin < 5
    rel_tol = 1.0e-6;
    if nargin < 4
        abs_tol = 1.0e-6;    
        if nargin < 3
            imax = 50;
        end
    end
end
[i,j] = size(X0);
if i<j
    X = X0.';
else 
    X = X0;
end
New_res_vec = zeros(imax,1);
New_abs_vec = zeros(imax,1);
New_rel_vec = zeros(imax,1);
Gmres_rel_vec = zeros(imax,1);
Gmres_flag_vec = zeros(imax,1);
test = true;
New_iter = 0;
u0 = X0(1:end-1);
u_theta_0 = dtheta_electro(u0,K,Nc);

while (test)
    unew = X(1:end-1);
    omega = X(end);
    phi_total =  Electro_time_stepper(tspan,unew,Ra,dt,Pr,Re,rr,K,Nc);
    phi = phi_total(:,end);
    gamma = rotation_symmetry(unew,omega,tspan,Nc,K);
    b = -[(phi-gamma) ; ...
          u_theta_0.'*(unew-u0)];
    Ax = @(v) matvec_prod_rel_equi(unew,phi,u_theta_0,omega,Ra,v);
    [delta,flag_gmres,relres,gmres_iter,resvec] = gmres(Ax,b,10,1e-3);
%     [delta,flag,relres,gmres_iter,resvec] = gmres(@(v) matvec_prod_parametization(unew,phi,c,Ra,v),b,10,1e-3);
%     [delta] = gmres(@(v) matvec_prod_parametization(unew,phi,c,Ra,v),b,30,1e-4);
    X = X + delta;
    New_res = norm(b,'inf');
    abs_err = norm(delta);
    if unew < 1e-6
        rel_err = abs_err;
    else
        rel_err = abs_err/norm(X);
    end
    root = X;
    New_iter = New_iter + 1;
    test =((New_res > abs_tol) & (abs_err > abs_tol) & (rel_err > rel_tol) ) & (New_iter < imax);
%     test =((New_res>1e-8) & (New_iter < imax));  
    disp(['parameter = ' num2str(Ra) '        ' ' Newton iteration = ' num2str(New_iter)])
    disp([' Newton residual = ' num2str(New_res),' abs error = ' num2str(abs_err), '    ' 'rel error = ' num2str(rel_err)])
    disp(['flag = ' int2str(flag_gmres),'    ' ' gmres_iter = ' int2str(gmres_iter),'    ' ' residual = ' num2str(relres)])
%     figure(New_iter)
%     semilogy(resvec),hold on
%     xlabel('GMRES Iteration','Interpreter', 'latex','Fontsize',16)
%     ylabel('Relative Residual','Interpreter', 'latex','Fontsize',16)
%     title('Convergence of GMRES residual','Interpreter', 'latex','Fontsize',16)
%     semilogy(New_iter ,New_res,'o',New_iter ,New_res,'*'),hold on 
%     xlabel('Newton iteration','Interpreter', 'latex','Fontsize',16)
%     ylabel('Newton residual','Interpreter', 'latex','Fontsize',16)
%     title('Convergence of Newton residual','Interpreter', 'latex','Fontsize',16)
%     pause
    New_res_vec(New_iter) = New_res;
    New_rel_vec(New_iter) = rel_err;
    New_abs_vec(New_iter) = abs_err;
    Gmres_rel_vec(New_iter) = relres;
    Gmres_flag_vec(New_iter) = flag_gmres;
end
tic
u_star = root(1:end-1);
angular_speed = root(end);
thetap = angular_speed*tspan(2);
omega_t = (thetap+2*pi/6)/tspan(2);
period = (2*pi/6)*(1/omega_t);
phi_period =  Electro_time_stepper([0, period],u_star,Ra,dt,Pr,Re,rr,K,Nc);
phi_period = phi_period(:,end);

[V,D,flag_eigs] = eigs(@(v) matvec_prod_rel_equi_eig(u_star,phi_period,Ra,period,v),4*(Nc+1)*(2*K+1),10);    
mu = diag(D);
lambda = log(diag(D))/tspan(2); 
disp(['eigs flag = ' , int2str(flag_eigs)])
disp([lambda mu])
eig_time = toc;
% figure()
% plot(real(mu),imag(mu),'bs'),grid on 
% title(['parameter = ' num2str(4/nu)])
% pause
%
% disp([log(diag(D))/tspan(2),diag(D)]);

      