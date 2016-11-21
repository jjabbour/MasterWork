init;
clf;
load('data_bifurcation_diagram.mat')
stable_periodic =[]; unstable_periodic = [];
for k = 1:size(periodic_result,2)
    if periodic_result(1,k) < 636
        stable_periodic = [stable_periodic, periodic_result(:,k)];
    else
        unstable_periodic = [unstable_periodic,periodic_result(:,k)];
    end
end

stable_parazero = 550:572;
unstable_parazero = 573:max(periodic_result(1,:));
stable_norm_zero = zeros(length(stable_parazero),1);
unstable_norm_zeros = zeros(length(unstable_parazero),1);

fig = figure(1);
plot(stable_periodic(1,:),stable_periodic(2,:),'o','MarkerFaceColor',[.1,.1,.1]),hold on
plot(unstable_periodic(1,:),unstable_periodic(2,:),'o')
plot(stable_parazero,stable_norm_zero,'-',unstable_parazero,unstable_norm_zeros,'b--','LineWidth',2)
plot([572.5,635.5],[0,mean([periodic_result(2,64),periodic_result(2,63)])],'rs','LineWidth',2)
%title('Bifurcation Diagram','Interpreter','latex','FontSize',16)
axis([550, 680, -.1 1.8])
xlabel('$\mathcal{R}$','Interpreter','latex','FontSize',20)
set(gca,'XTick',[(min(stable_parazero)),572,590,600,636, max(periodic_result(1,:))],'FontSize',14)
set(gca,'XTickLabel',{int2str(min(stable_parazero)),'572','590','600','636',int2str(max(periodic_result(1,:)))})
ylabel('$\|\phi\|_{2}$ ', 'Interpreter','latex','FontSize',20)
h_leg = legend('stable limit cycles','unstable limit cycles',...
        'stable steady states','unstable steady states','bifurcation point','Location','NorthWest');
set(h_leg,'FontSize',16,'Location','NorthWest');
% save the plot as a pdf file
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition',[0 0 1 1]);
% print(fig,'-dpdf','Bifurcation_Diagram')   
print(fig,'-dpdf','..\Pictures\Bifurcation_Diagram.pdf') 