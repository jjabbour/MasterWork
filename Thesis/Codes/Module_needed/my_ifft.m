function U = my_ifft(um,N)
% If I have K+1 Fourier coefficients for the zero mode and first K modes of
% a real-valued function, I want to use these to construct a real-valued
% function, U.
%
%
% if u is on 64 meshpoints then fft(u) will have the zero mode in the first
% entry, will have modes 1 to 31 in the 2 to 32 entries, will have scruff
% in the 33rd entry, and will have modes -31 to -1 in the 34 to 64th
% entries.  NOTE: if I want fft and ifft to be exact inverses of one
% another, I cannot throw away that scruff.  It might not be at the level
% of roundoff error...  
%
% If u is on 63 mesh points then fft(u) will have the zero mode in the first
% entry, will have modes 1 to 31 in the 2 to 32 entries, and will have 
% modes -31 to -1 in the 33 to 63th entries.
K = ceil(N/2-1);

um=transpose(um);

Um(1,:)=um(1,:);
if mod(N,2)==0
    % N is even
    Um(2:K+2,:)=um(2:K+2,:);
    Um(K+3:N,:)=conj(um(K+1:-1:2,:));
else
    % N is odd
    Um(2:K+1,:)=um(2:K+1,:);
    Um(K+2:N,:)=conj(um(K+1:-1:2,:));   
end
U = transpose(ifft(Um));
% 
% N = 7
% K = ceil(N/2-1)
% if mod(N,2)==0
%     % N is even
%     2:K+2
%     K+3:N
%     K+1:-1:2
% else
%     % N is odd
%     2:K+1
%     K+2:N
%     K+1:-1:2  
% end