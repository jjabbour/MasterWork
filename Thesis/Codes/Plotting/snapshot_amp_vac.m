init
load('time_stepping1em1_Ra_628  632  640  644.mat')
close all;
U_Ra_640 = U_per(:,:,3);
a = 10000; step = 750;
plot_vec = a:step:a+3*step;
path_name = 'C:\Users\100380000\Desktop\Master_Thesis\Thesis_paper\figures\Pictures\';
% path_name = '..\Pictures\'
for k = 1:4
    [a,b,c,d] = vec_2_mat(U_Ra_640(:,plot_vec(k)),Nc,K);
%     subplot(2,2,k)
    fig1 = figure(k); 
    polar3d(d,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
    set(gca,'FontSize',16)
    set(fig1,'PaperOrientation','landscape');
    set(fig1,'PaperUnits','normalized');
    set(fig1,'PaperPosition',[0 0 1 1]);
    print(fig1,'-dpdf',[path_name,'snapshot_amp_vac', int2str(k),'.pdf'])  

end
    
% for ell = 1:size(U_per,2)
%     U_plot = U_per(:,:,ell);
%     for k = 1:100:size(U_plot,2)
%         [a,b,c,d] = vec_2_mat(U_plot(:,k),Nc,K);
%         polar3d(d,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
%         title(['Ra = ' int2str(Ra(ell))],'interpreter','latex','FontSize',20)
%         pause(.01)
%     end
% end