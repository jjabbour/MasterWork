init
% load('time_stepping1em1_Ra_628  632  640  644.mat')

for ell = 1:size(U_per,2)
    U_plot = U_per(:,:,ell);
    figure(1); clf; pause
    for k = 1:100:size(U_plot,2)
        [a,b,c,d] = vec_2_mat(U_plot(:,k),Nc,K);
        polar3d(d,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');
        title(['Ra = ' int2str(Ra(ell))],'interpreter','latex','FontSize',20)
        pause(.01)
    end
end
