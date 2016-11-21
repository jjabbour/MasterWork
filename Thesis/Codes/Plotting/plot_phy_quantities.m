init
load('tracing_periodic_solution_date735924');
col = 53;
% sort([X(end,:) ;abs(mu)],'descend');
U_star = X(1:end-2,col);
par_star = X(end,col);
[psi_2,q,omega,phi] = vec_2_mat(U_star,Nc,K);

fig1 = figure(1);
polar3d(psi_2,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
set(gca,'FontSize',16)
title('(a)                                       ','FontSize',20)
%title('Perturbation of the 2D electric potential difference $\psi_2$','Interpreter','Latex','FontSize',24)
fig2 = figure(2);
polar3d(q,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
set(gca,'FontSize',16)
title('(b)                                       ','FontSize',20)
%title('Perturbation of the surfce charge density $q$ ','Interpreter','Latex','FontSize',24)
fig3 = figure(3);
polar3d(omega,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
set(gca,'FontSize',16)
%title('Perturbation of the vorticity $\omega$','Interpreter','Latex','FontSize',24)
title('(c)                                       ','FontSize',20)
fig4 = figure(4);
polar3d(phi,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
set(gca,'FontSize',16)
%title('Perturbation of the streamfunction $\phi$ ','Interpreter','Latex','FontSize',24)
title('(d)                                       ','FontSize',20)
path_name = 'C:\Users\100380000\Desktop\Master_Thesis\Thesis_paper\figures\Pictures\';
% path_name = '..\Pictures\';

set(fig1,'PaperOrientation','landscape');
set(fig1,'PaperUnits','normalized');
set(fig1,'PaperPosition',[0 0 1 1]);
print(fig1,'-dpdf',[path_name 'pot_ele.pdf']);

set(fig2,'PaperOrientation','landscape');
set(fig2,'PaperUnits','normalized');
set(fig2,'PaperPosition',[0 0 1 1]);
print(fig2,'-dpdf',[path_name 'charge_dens.pdf']);

set(fig3,'PaperOrientation','landscape');
set(fig3,'PaperUnits','normalized');
set(fig3,'PaperPosition',[0 0 1 1]);
print(fig3,'-dpdf',[path_name 'vorticity.pdf']);

set(fig4,'PaperOrientation','landscape');
set(fig4,'PaperUnits','normalized');
set(fig4,'PaperPosition',[0 0 1 1]);
print(fig4,'-dpdf',[path_name 'stream_fun.pdf']);
