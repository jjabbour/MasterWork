function um = my_fft(U,N)
% If I have a real-valued function at N points, I want to extract the
% coefficients of the zero mode and the K positive modes.
K = ceil(N/2 - 1);

Um = fft(transpose(U));
size(Um);
if mod(N,2)==0
    % N is even, need to keep the scruff mode
    um = Um(1:K+2,:);
else
    % N is odd, there's no scruff mode
    um = Um(1:K+1,:);
end
um = transpose(um);


