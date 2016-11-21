init
% load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_5Ra_637to647')
% load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_5Ra_636to639')
% load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_10Ra_636to639')
load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_30Ra_635to646')
% load('direct_simulation_of_rotating_wave_with_pert_after_bifu_for_t_30Ra_638to638.mat')
[m,n,p] = size(sim_output);
% t = linspace(tf-5,tf,n); % when you store the last 10000 output
t = linspace(0,tf,n); 


for k = 1:p
%     tic
    hump = [];
    x_temp   = sim_output(5115,:,k);
    vec_diff = sign(diff(x_temp));
    for ell = 1: length(vec_diff)-1
        test = vec_diff(ell)>vec_diff(ell+1); 
        if test
            hump = [hump,x_temp(ell+1)];
        end
    end
    uu = 8;
    figure(1)
    plot(t,x_temp)
    for ii = 1:floor(length(hump)/uu)
        max_amp(ii) = max(hump(uu*(ii-1)+1:uu*ii));
        min_amp(ii) = min(hump(uu*(ii-1)+1:uu*ii));
    end
    amp_vac = 1/2*(max_amp - min_amp);
    figure()
    plot(amp_vac,'*'),axis([0 length(amp_vac)  min(amp_vac) max(amp_vac)])
end

% figure(p+1)
% plot(Ra_v,amp_vac,'r-o'),
% % axis([min(Ra_v) max(Ra_v) min(amp_vac)-1e-1 max(amp_vac)+1e-1])
% xlabel('$\mathrm{Ra}$','Interpreter','Latex','FontSize',14)
% ylabel('Amplitude difference','FontSize',14)

% figure(p+2)
% plot(Ra_v,norm_phi,'bo')
% xlabel('$\mathrm{Ra}$','Interpreter','Latex','FontSize',14)
% ylabel('$\mathrm{norm}\{\phi(r,\theta,t)\}$','Interpreter','Latex','FontSize',14)

