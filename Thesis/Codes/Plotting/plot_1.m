close all
init
load('time_stepping_the_two_solution_bef_aft_second_bif');
color = ['m' 'b' 'g' 'y' 'k'];
x1_star_Ra_632 = Up_star_Ra_632(5115,1:end);
x2_star_Ra_632 = Up_star_Ra_632(6165,1:end);

x1_star_Ra_640 = Up_star_Ra_640(5115,1:end);
x2_star_Ra_640 = Up_star_Ra_640(6165,1:end);
x3_star_Ra_640 = Up_star_Ra_640(6345,1:end);

% plot(x1,x2,'r'), hold on 
% plot(x1_star,x2_star,'r'), hold on 
%plot(y1_star,y2_star,'r'), hold on
fig1 = figure(1);
plot(x1_star_Ra_632,x2_star_Ra_632); hold on 
xlabel('$x_1$','Interpreter','latex','FontSize',14)
ylabel('$x_2$','Interpreter','latex','FontSize',14)
title(['Ra = ' int2str(Ra_v(1))],'Interpreter','latex','FontSize',20)
fig2 = figure(2);
plot(x1_star_Ra_640,x2_star_Ra_640);hold on 
xlabel('$x_1$','Interpreter','latex','FontSize',14)
ylabel('$x_2$','Interpreter','latex','FontSize',14)
title(['Ra = ' int2str(Ra_v(2))],'Interpreter','latex','FontSize',20)
for k = 1:length(color)
    figure(1)
    plot(x1_star_Ra_632(end-k*(2000):end-(k-1)*2000),x2_star_Ra_632(end-k*(2000):end-(k-1)*2000),...
        color(k),'LineWidth',2);
    figure(2)
    plot(x1_star_Ra_640(end-k*(2000):end-(k-1)*2000),x2_star_Ra_640(end-k*(2000):end-(k-1)*2000),...
        color(k),'LineWidth',2);
end
% plot time series 
dt = 5e-4;
tf = size(Up_star_Ra_640,2)*5e-4;
t = 0:dt:tf; t = t(1:end-1);
u1_star_Ra_632 = Up_star_Ra_632(5115,:);
u1_star_Ra_640 = Up_star_Ra_640(5115,:);
mn = length(t); a =8000;st = mn-a; 
fig3 = figure(3);
% plot(t,x1_star_Ra_632(end-length(t)+1:end))
plot(t(st:mn),u1_star_Ra_632(st:mn)), axis([t(st) t(mn) -.12 .12])
fig4 = figure(4);
% plot(t,x1_star_Ra_640(end-length(t)+1:end))
plot(t(st:mn),u1_star_Ra_640(st:mn)),axis([t(st) t(mn) -.12 .12])


set(fig1,'PaperOrientation','landscape');
set(fig1,'PaperUnits','normalized');
set(fig1,'PaperPosition',[0 0 1 1]);
print(fig1,'-dpdf','..\Pictures\Sim_Ra_632.pdf')

set(fig2,'PaperOrientation','landscape');
set(fig2,'PaperUnits','normalized');
set(fig2,'PaperPosition',[0 0 1 1]);
print(fig2,'-dpdf','..\Pictures\Sim_Ra_640.pdf')

set(fig3,'PaperOrientation','landscape');
set(fig3,'PaperUnits','normalized');
set(fig3,'PaperPosition',[0 0 1 1]);
print(fig3,'-dpdf','..\Pictures\timeserie_632.pdf')

set(fig4,'PaperOrientation','landscape');
set(fig4,'PaperUnits','normalized');
set(fig4,'PaperPosition',[0 0 1 1]);
print(fig4,'-dpdf','..\Pictures\timeserie_640.pdf')