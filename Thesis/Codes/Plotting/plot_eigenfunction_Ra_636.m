init
load('tracing_periodic_solution_date735924')
xk = 57 ; 
% rename the data of the computation 
X_star = X(:,xk);
Eigen_vector = V(:,:,xk); 
lambda = mu(:,xk);

U_star = X_star(1:end-2);
omega_star = X_star(end-1);
Ra  = X_star(end);
ell = 1;
for k = 1:10
    if abs(lambda(k)) > 1
        A(:,ell)  = [lambda(k); Eigen_vector(:,k)];
        ell = ell + 1;
    end
end

Eigen_function = A(2:end,:);

for h = 1:size(A,2)
    fig1 = figure(h); 
    [a,b,c,d] = vec_2_mat(real(Eigen_function(:,h)),Nc,K);
    polar3d(d,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
    title(['$\lambda$ = ' num2str(A(1,h)) , ' and  $|\lambda|$ = ' , num2str(abs(A(1,h)))],...
        'Interpreter', 'Latex','FontSize',20 )
    
    pic_name = ['Eigen_function_', int2str(h) '.pdf'];
    set(fig1,'PaperOrientation','landscape');
    set(fig1,'PaperUnits','normalized');
    set(fig1,'PaperPosition',[0 0 1 1]);
    print(fig1,'-dpdf',['..\Pictures\', pic_name])
end
 


