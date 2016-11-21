function [psi2,psi2m,q,qm,w,wm,phi,phim] = random_initial_data(Tq_vs_psi2,D,D2,riprime,theta,Nc,K,InitBCs_amp)
% Random Noise as initial conditions: 
% specify the initial conditions as random noises for stream and vorticity function 
% All Phi = 0, w=0, but Psi2 = random noise * InitBCs_amp for all Fourier modes

% to populate all frequencies w/ noise:
psi2 = zeros(Nc+1,2*K+1);
for i=2:Nc
    for j=0:8
        psi2(i,:) = psi2(i,:)+rand*cos(j*(theta+2*pi*rand));
    end
end
% specify the amplitude...
psi2 = psi2/max(max(psi2));
psi2 = psi2*InitBCs_amp;
% now multiply by something to make it smooth in the radial direction...
% and to make it zero at the inner and outer boundaries...
window = (riprime(1)-riprime).^3.*(riprime-riprime(Nc+1)).^3;
window = window/max(window);
psi2 = diag(window)*psi2;
% to get initial data for charge density
psi2m = my_fft(psi2,2*K+1);
for m=0:1:K              %  Assign \hat{q}_m in accordance with the values of  \hat{psi2}_m.
    qm(:, m+1)=Tq_vs_psi2(:,:, m+1)*psi2m(:, m+1);
end
q = my_ifft(qm,2*K+1);

phi = zeros(Nc+1,2*K+1);
phim = my_fft(phi,2*K+1);
w = zeros(Nc+1,2*K+1);
wm = my_fft(w,2*K+1);

% for testing purposes...
% phim = rand(size(phim));
% to populate all frequencies w/ noise:
phi = zeros(Nc+1,2*K+1);
for i=2:Nc
    for j=0:8
        phi(i,:) = phi(i,:)+rand*cos(j*(theta+2*pi*rand));
    end
end
% specify the amplitude...
phi = phi/max(max(phi));
phi = phi*InitBCs_amp;
% now multiply by something to make it smooth in the radial direction...
% and to make it zero at the inner and outer boundaries...
phi = diag(window)*phi;
phim = my_fft(phi,2*K+1);

w = zeros(Nc+1,2*K+1);
for m=0:1:K
    % this is the Laplacian
    Lop = 4*D2 + 2*diag(1./riprime)*D - m^2*diag(1./(riprime.^2));
    wm(:,m+1) = -Lop*phim(:,m+1);
end
w = my_ifft(wm,2*K+1);

% for i=2:Nc
%     for j=0:8
%         w(i,:) = w(i,:)+rand*cos(j*(theta+2*pi*rand));
%     end
% end
% % specify the amplitude...
% w = w/max(max(w));
% w = w*InitBCs_amp;
% % now multiply by something to make it smooth in the radial direction...
% % and to make it zero at the inner and outer boundaries...
% w = diag(window)*w;
% wm = my_fft(w,2*K+1);