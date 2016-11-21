function [psi2m_last,qm_last,wm_last,phim_last,Um_init] = reshape_IC_time_stepper(unew,Nc,K)
% [psi2,q,w,phi] 2D potential , charge, vorticity, stream line, each (2K+1)x(Nc+1)
% Input 
% unew  4(2K+1)(Nc+1) x  1   
% Nc    1 x 1                   Highest order of Chebychev polynomial
% K     1 x 1                   Highest wave number of Fourrier bases function 

% Output:
% Um_init   4(K+1)(Nc+1) x 1    Spectral transform of unew to be initial conditiion for the time stepper implemented 

m = length(unew);
q = 4*(2*K+1)*(Nc+1);
if m ~= q
    unew = unew(1:q,1);
end
    
[psi2_last,q_last,w_last,phi_last] = vec_2_mat(unew,Nc,K);
psi2m_last = my_fft(psi2_last,2*K+1);   U1m = reshape(psi2m_last,(K+1)*(Nc+1),1);
qm_last = my_fft(q_last,2*K+1);         U2m = reshape(qm_last,(K+1)*(Nc+1),1);
wm_last = my_fft(w_last,2*K+1);         U3m = reshape(wm_last,(K+1)*(Nc+1),1);
phim_last = my_fft(phi_last,2*K+1);     U4m = reshape(phim_last,(K+1)*(Nc+1),1);
Um_init = [U1m;U2m;U3m;U4m];

end