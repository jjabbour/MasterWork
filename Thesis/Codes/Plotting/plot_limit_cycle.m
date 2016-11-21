init
load('result_735895_in_increasing _order.mat')
U = X(1:end-2,4);
phase_speed = X(end-1,4);
Ra  = X(end,4);
period_p = period(4);

load('tracing_periodic_solution_date735924.mat')
points = 6:10:26;

U = [U X(1:end-2,points)];
phase_speed = [phase_speed X(end-1,points)];
Ra  = [Ra X(end,points)];
period_p = [period_p period(points)];

for k = 1:length(points)+1
    [Up,Uf] = Electro_periodic_check([0 period(k)],U(:,k),Ra(k));
    x1 = Up(5168,:);
    x2 = Up(6208,:);
    [a,b,c,d] = vec_2_mat(Up(:,end),Nc,K);
    fig1 = figure(1);
    subplot(2,2,k)
    polar3d(d,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
    title({'Perturbation of the streamfunction $\phi(r,\theta)$ ';['from the base state at Ra = ' int2str(Ra(k))]},'Interpreter','Latex','FontSize',14)
    title(['Ra = ' int2str(Ra(k))],'Interpreter','Latex','FontSize',14)
    fig2 = figure(2);
    subplot(2,2,k)
    plot(x1,x2)
    xlabel('$u(r_1,\theta_1)$','Interpreter','Latex','Fontsize', 12)
    ylabel('$u(r_2,\theta_2)$','Interpreter','Latex','Fontsize', 12)
%     title(['Limit cycle  for Ra = ', int2str(Ra(k))],'Interpreter','Latex','Fontsize', 14)
    title(['Ra = ', int2str(Ra(k))],'Interpreter','Latex','Fontsize', 14)
    legend(['period = ' num2str(period(k),'%.6f')],'Location','NorthWest')   
end

% % saving figure 1
% set(fig1,'PaperOrientation','landscape');
% set(fig1,'PaperUnits','normalized');
% set(fig1,'PaperPosition',[0 0 1 1]);
% print(fig1,'-dpdf','..\Pictures\perturbations_periodic')
% % print(fig1,'-dpdf','C:\Users\100380000\Desktop\Master_Thesis_codes\Pictures\Perturbation_periodic.pdf')  
% 
% % saving figure(2) 
% set(fig2,'PaperOrientation','landscape');
% set(fig2,'PaperUnits','normalized');
% set(fig2,'PaperPosition',[0 0 1 1]);
% % print(fig2,'-dpdf','limit_cycle')  
% print(fig2,'-dpdf','..\Pictures\limit_cycle.pdf')      
    