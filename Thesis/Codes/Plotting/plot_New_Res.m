init
load('tracing_periodic_solution_date735924')
clf
points = 2:5:27;
bif_param = X(end,points);
% legend_string = zeros(length(points),19);
fig = figure(1);
for k = 1:length(points)
    legend_string(k,:) = ['$\mathcal{R}$ = ', int2str(bif_param(k)) ];
end
% figure()
semilogy(NEW_res(:,points),'-o')
xlabel('$k$','Interpreter','Latex','FontSize',20)
ylabel('Residual' ,'Interpreter','Latex','FontSize',20)
%title('Convergence of the  nonlinear solver' ,'Interpreter','Latex','FontSize',20)
h_leg = legend(legend_string);
set(h_leg,'FontSize',20,'Interpreter','Latex','Location','NorthEast')
set(gca,'FontSize',14)
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition',[0 0 1 1]);
print(fig,'-dpdf','..\Pictures\Newton_residual.pdf')  

