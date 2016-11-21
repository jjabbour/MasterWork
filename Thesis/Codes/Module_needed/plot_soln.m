
subplot(2,2,1)
polar3d(psi2_phy(:,:,ii+1),0,2*pi,rr,1,1,'contour');
title('Psi2')

subplot(2,2,2)
polar3d(q_phy(:,:,ii+1),0,2*pi,rr,1,1,'contour');
title('Q');

subplot(2,2,3)
polar3d(Phi_phy(:,:,ii+1),0,2*pi,rr,1,1,'contour');
title('Phi')

subplot(2,2,4)
polar3d(W_phy,0,2*pi,rr,1,1,'contour');
title('W')


% figure(3)
% polar3d(W_phy,0,2*pi,min(riprime),max(riprime),1,'contour');
% theta_plot = theta;
% theta_plot(2*K+2) = 2*pi;
% [rg,thetag]=meshgrid(riprime,theta_plot);
% X=rg.*cos(thetag);
% Y=rg.*sin(thetag);
% W_phy(:,2*K+2) = W_phy(:,1);
% contourf(X,Y,transpose(W_phy));
% colorbar
% axis equal
% % Clim sets the lower and upper bounds, it lets matlab pick colours for the
% % contour plot.
% set(gca, 'Clim',[-max(riprime),max(riprime)])