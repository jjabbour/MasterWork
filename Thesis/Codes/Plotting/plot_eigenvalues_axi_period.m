init
points = 1:30;

n = 10000000;
k = 1:n; 
xs = cos(2*pi*k/n);
ys = sin(2*pi*k/n);

% load('tracing_zero_solution_branch735917.mat')
lambda = mu(:,points);
Ra_v   = X(end,points);

fig = figure(1);

a = 1.55;
axis([-a a -a a])
axis('square')
grid('on') 
hold on 
plot(real(lambda(:,1)),imag(lambda(:,1)),'bs')
plot(real(lambda(:,23)),imag(lambda(:,23)),'ko','MarkerFaceColor',[0,0,0]);
plot(real(lambda),imag(lambda),'r.')
plot(xs,ys,'k'),
%plot(real(lambda(:,38)),imag(lambda(:,38)),'blacko');

xlabel('Re $\,\lambda$','Interpreter','Latex','FontSize',20)
ylabel('Im $\,\lambda$','Interpreter','Latex','FontSize',20)
% title({'The evolution of the largest 10 eigenvalues of the linearization of the map $\Phi$ ' ;...
%     [' about the axisymetric solution for Rayleigh Ra between ' int2str(min(Ra_v)),' and Ra = ' int2str(max(Ra_v))]},...
%     'Interpreter','Latex','FontSize',20); 
% title('Eigenvalues of $\mathrm{D_u}\Phi_{t}$',...
%         'Interpreter','Latex','FontSize',20); 
set(gca,'FontSize',14)
h_leg = legend(['$\mathcal{R}$ = ',int2str(Ra_v(1))],...
        ['$\mathcal{R}$ = ',int2str(Ra_v(23))],'$\lambda(\mathcal{R})$');
set(h_leg,'FontSize',16,'Interpreter','Latex','Location','NorthEast')    

% legend(['Ra = ',int2str(Ra_v(1))], ['Ra = ',int2str(Ra_v(23))],'evolution of eigenvalues') 

% % save the plot as a pdf file
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition',[0 0 1 1]);
% print (fig,'-dpdf','Eigenvalue_axi_perio')   
print(fig,'-dpdf','..\Pictures\Eigenvalue_axi_perio.pdf') 
