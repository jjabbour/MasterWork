init
load('plot_axisymetrix_to_peridic.mat')

K = 32; Nc =24; rr =.56;

[a0,b0,c0,d0] = vec_2_mat(Upco(:,3),Nc,K);
[af,bf,cf,df] = vec_2_mat(Uf1,Nc,K);
fig = figure(1);
subplot(1,2,1)
polar3d(df,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
% title({'Axisymmetric Solution'},'Interpreter','latex','FontSize',20,'FontWeight','bold')
title('(a)                               ','FontSize',20)
subplot(1,2,2)
polar3d(d0,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
% title({'Rotating Wave '},'Interpreter','latex','FontSize',20,'FontWeight','bold')
title('(b)                               ','FontSize',20)
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition',[0 0 1 1]);
print(fig,'-dpdf','..\Pictures\axi_rot_plot.pdf')


