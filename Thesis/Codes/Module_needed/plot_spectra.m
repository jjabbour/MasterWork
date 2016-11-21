    Psi2_spec = log10(abs(psi2m_now)+10^(-16));
    Q_spec = log10(abs(qm_now)+10^(-16));
    Phi_spec = log10(abs(phim_now)+10^(-16));
    W_spec = log10(abs(wm_now)+10^(-16));
    
    max(max(abs(phim_now)));

    figure(1)
    subplot(2,2,1)
    plot(Psi2_spec')
    title('Psi2 spec')
    axis([0,K+1,-17,2])
    subplot(2,2,2)
    plot(Q_spec')
    title('Q spec');
    axis([0,K+1,-17,2])
    subplot(2,2,3)
    plot(Phi_spec')
    title('Phi spec')
    axis([0,K+1,-17,2])
    subplot(2,2,4)
    plot(W_spec')
    title('W spec')
    axis([0,K+1,-17,2])
    figure(1)
    % pause(1)