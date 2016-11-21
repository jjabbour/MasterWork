clf;
init
points = 30:70;

n = 10000000;
k = 1:n; 
xs = cos(2*pi*k/n);
ys = sin(2*pi*k/n);
% load('tracing_periodic_solution_date735924.mat')


lambda = mu(:,points);
Ra_v   = X(end,points);

fig = figure(1);
a = 1.55;
axis([-a a -a a])
axis('square')
grid on;
hold on; 
plot(real(lambda(:,1)),imag(lambda(:,1)),'bs ')
plot(real(mu(:,57)),imag(mu(:,57)),'ko','MarkerFaceColor',[0,0,0]);
plot(real(lambda),imag(lambda),'r.')
plot(xs,ys,'k'),

xlabel('Re$\,\lambda$','Interpreter','Latex','FontSize',20)
ylabel('Im$\,\lambda$','Interpreter','Latex','FontSize',20)
% title({'The evolution of the largest 10 eigenvalues of the linearization of the map $\Phi$ ' ; [' about a periodic '...
%     ' solution for Rayleigh Ra between ' int2str(min(Ra_v)),' and Ra = ' int2str(max(Ra_v))]},...
%     'Interpreter','Latex','FontSize',20); 
% title('Eigenvalues of $\mathrm{D_u}\Phi_{\tau}$ ',...
%       'Interpreter','Latex','FontSize',20); 
set(gca,'FontSize',14)
h_leg = legend(['$\mathcal{R}$ = ',int2str(Ra_v(1))],...
        ['$\mathcal{R}$ = ',int2str(X(end,57))],'$\lambda(\mathcal{R})$');
set(h_leg,'FontSize',16,'Interpreter','Latex','Location','NorthEast')    

% save the plot as a pdf file
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition',[0 0 1 1]);
print(fig,'-dpdf','..\Pictures\Eigenvalue_period_torus.pdf') 
