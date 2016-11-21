init
% load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_5Ra_637to647')
% load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_5Ra_636to639')
% load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_10Ra_636to639')
% load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_30Ra_635to646')
% close all;
[m,n,p] = size(sim_output);
t = linspace(tf-5,tf,n);
% A = permute(sim_output,[2,1,3]);
max_amp = zeros(1,p);
min_amp = zeros(1,p);

for k = 1:p
    hump = [];
    x_temp   = sim_output(5115,:,k);
    vec_diff = sign(diff(x_temp));
    for ell = 1: length(vec_diff)-1
        test = vec_diff(ell)>vec_diff(ell+1); 
        if test
            hump = [hump,x_temp(ell+1)];
        end
    end
    max_amp(k) = max(hump(end-10:end));
    min_amp(k) = min(hump(end-10:end));
end

% for ell = 2:4:length(Ra_v)
%     figure()
%     plot(t,sim_output(5115,:,ell),'b-'), hold on
%     plot(t,max_amp(ell),'r-'),axis([tf-5 tf -max_amp(ell)-.5e-1 max_amp(ell)+.5e-1])
%     xlabel('$t$','Interpreter','Latex','FontSize',14)
%     ylabel('$\phi(r_{i},\theta_{i},t)$','Interpreter','Latex','FontSize',14)
%     title(['Ra = ',int2str(Ra_v(ell))],'FontSize',14)
%     legend('\phi(r_{i},\theta_{i},t)','max','Location','NorthWest')
% end

amp_vac = 1/2*(max_amp - min_amp);
fig1 = figure(1);
% pause
plot(Ra_v,amp_vac,'ro','MarkerFaceColor',[1,0,0]),
axis([min(Ra_v) max(Ra_v) min(amp_vac)-1e-1 max(amp_vac)+1e-1])
xlabel('$\mathcal{R}$','FontSize',14,'Interpreter','Latex')
ylabel('The amplitude of the vacillating wave','FontSize',14,'Interpreter','Latex')

set(fig1,'PaperOrientation','landscape');
set(fig1,'PaperUnits','normalized');
set(fig1,'PaperPosition',[0 0 1 1]);
print(fig1,'-dpdf','..\Pictures\Amplitude_vacillating.pdf')  

