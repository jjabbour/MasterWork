init;
load('time_stepping_the_two_solution_bef_aft_second_bif');

x1_star_Ra_632 = Up_star_Ra_632(5115,1:end);
x2_star_Ra_632 = Up_star_Ra_632(6165,1:end);

x1_star_Ra_640 = Up_star_Ra_640(5115,1:end);
x2_star_Ra_640 = Up_star_Ra_640(6165,1:end);
x3_star_Ra_640 = Up_star_Ra_640(6245,1:end);
dt = 5e-4;
tf = size(Up_star_Ra_640,2)*dt;
t = 0:dt:tf; t = t(1:end-1);
a = 8000; mn = length(x1_star_Ra_632); st = mn-a; k = 1;
t_plot = t(st:mn); 
u1_Ra_632 = x1_star_Ra_632(st:mn);
u2_Ra_632 = x2_star_Ra_632(st:mn);
u1_Ra_640 = x1_star_Ra_640(st:mn);

fig1 = figure(1);
fig2 = figure(2);
fig3 = figure(3);
fig4 = figure(4);

figure(1)
%subplot(2,2,1)
plot(t_plot,u1_Ra_632), axis([t(st) t(mn) -.12 .12])
xlabel('t','Interpreter','Latex','FontSize',14)
ylabel('$u_{1}(r,\theta,t)$','Interpreter','Latex','FontSize',14)
%title('(a)                                           ','FontSize',18)
%subplot(2,2,3)
figure(2)
plot(u1_Ra_632,u2_Ra_632);axis([-.12 .12 -.12 .12]);
xlabel('$u_{1}(r,\theta,t)$','Interpreter','Latex','FontSize',14)
ylabel('$u_{2}(r,\theta,t)$','Interpreter','Latex','FontSize',14)
%title('(b)                                           ','FontSize',18)
%subplot(2,2,2)
figure(3)
plot(t_plot,u1_Ra_640),axis([t(st) t(mn) -.12 .12])
xlabel('t','Interpreter','Latex','FontSize',14)
ylabel('$u_{1}(r,\theta,t)$','Interpreter','Latex','FontSize',14)
% title('(c)                                           ','FontSize',18)
% subplot(2,2,4)
figure(4)
x1 = x1_star_Ra_640(end-k*a:end-(k-1)*a); 
x2 = x2_star_Ra_640(end-k*a:end-(k-1)*a);
plot(x1,x2);axis([-.12 .12 -.12 .12]);
xlabel('$u_{1}(r,\theta,t)$','Interpreter','Latex','FontSize',14)
ylabel('$u_{2}(r,\theta,t)$','Interpreter','Latex','FontSize',14)
% title('(d)                                          ','FontSize',18)

% figure()
% plot3(x1_star_Ra_640,x2_star_Ra_640,x3_star_Ra_640)
set(fig1,'PaperOrientation','landscape');
set(fig1,'PaperUnits','normalized');
set(fig1,'PaperPosition',[0 0 1 1]);
print(fig1,'-dpdf','..\Pictures\time_phase_Ra_632_640_1.pdf')

set(fig2,'PaperOrientation','landscape');
set(fig2,'PaperUnits','normalized');
set(fig2,'PaperPosition',[0 0 1 1]);
print(fig2,'-dpdf','..\Pictures\time_phase_Ra_632_640_2.pdf')

set(fig3,'PaperOrientation','landscape');
set(fig3,'PaperUnits','normalized');
set(fig3,'PaperPosition',[0 0 1 1]);
print(fig3,'-dpdf','..\Pictures\time_phase_Ra_632_640_3.pdf')

set(fig4,'PaperOrientation','landscape');
set(fig4,'PaperUnits','normalized');
set(fig4,'PaperPosition',[0 0 1 1]);
print(fig4,'-dpdf','..\Pictures\time_phase_Ra_632_640_4.pdf')