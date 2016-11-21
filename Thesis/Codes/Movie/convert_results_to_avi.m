init
%load('time_stepping1em1_Ra_628  632  640  644.mat')

[r_size,c_size,d_size] = size(U_per);
vec = 10:10:c_size;
nFrames = length(vec);


% Preallocate movie structure.
mov(1:nFrames) = struct('cdata', [],...
                        'colormap', []);

% Create movie.
% [a,b,c,d] = vec_2_mat(Up_star_Ra_640(:,30000),Nc,K);
for ell =  1:d_size
    [a,b,c,d] = vec_2_mat(U_per(:,10,ell),Nc,K);
    polar3d(d,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour');axis tight
    set(gca,'nextplot','replacechildren');
    for k = 1:nFrames 
        [a,b,c,d] = vec_2_mat(U_per(:,vec(k),ell),Nc,K);
        polar3d(d,0,2*pi,rr/(1-rr),1/(1-rr),1,'contour'); colormap
        title(['$\mathcal{R}$ = ',int2str(Ra(ell))],'Interpreter','latex','FontSize',20)
        pause%(.001)
        mov(k) = getframe(gcf);
    end
    % Create AVI file.
    movie2avi(mov, ['Ra' int2str(Ra(ell)) '.avi'], 'compression', 'None');
end

