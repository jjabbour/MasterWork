function   [wm1,phim1,wm3,phim3,Mi] = influence_matrix(D,D2,riprime,K,dt,Pr,It)


% ==========================================================
%  The Influence Matrix Method
%               time-independent problem for wm1 and phim1, wm3 and phim3
% ===========    Step 1. \bar{P}_1 and \bar{P}_2 problems  ============
[a,b] = size(D);
Nc=a-1
for m=0:1:K % Loop for each Fourier mode
    Lop = 4*D2 + 2*diag(1./riprime)*D - m^2*diag(1./(riprime.^2)); % Laplace operator matrix (checked)
    % now we want to solve
    %             (I-delta L) w = 0
    %                    L  phi = -w
    % subject to the boundary conditions phi(\pm 1) = 0 and w(-1) = 1, w(1) = 0
    %
    % first solve for w
    A_full = (It - 2/3*dt*Pr*Lop(:,:));
    A_trunc = A_full(2:Nc,2:Nc);
    wm1(1,m+1)=0;
    wm1(Nc+1,m+1)=1;
    B = 0;
    RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[wm1(1,m+1);wm1(Nc+1,m+1)];
    wm1(2:Nc,m+1)=A_trunc\RHS;
    % max(abs(A_full(2:Nc,:)*wm1(:,m+1)-B))
    %
    % now solve for phi
    A_full = Lop;
    A_trunc = A_full(2:Nc,2:Nc);
    phim1(1,m+1)=0;
    phim1(Nc+1,m+1)=0;
    B = -wm1(2:Nc,m+1);
    RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[phim1(1,m+1);phim1(Nc+1,m+1)];
    phim1(2:Nc,m+1)=A_trunc\RHS;
    %max(abs(A_full(2:Nc,:)*phim1(:,m+1)-B))
    
    %% Peichun's computation
    % wm1_PT(1:Nc+1, m+1) = ( It - 2/3*dt*Pr*[ zeros(1,Nc+1);Lop(2:Nc, :); zeros(1,Nc+1)] )\ [zeros(Nc,1); 1];
    % phim1_PT(2:Nc, m+1) = -  Lop(2:Nc,2:Nc)\ wm1_PT(2:Nc, m+1);
    % phim1_PT(1, m+1) = 0;
    % phim1_PT(Nc+1, m+1) = 0;
    % compare her solution to my solution:
    % max(abs(wm1(:,m+1)-wm1_PT(:,m+1)));
    %max(abs(phim1(:,m+1)-phim1_PT(:,m+1)))
    
    % now we want to solve
    %             (I-delta L) w = 0
    %                    L  phi = -w
    % subject to the boundary conditions phi(\pm 1) = 0 and w(-1) = 0, w(1) = 1
    %
    % first solve for w
    A_full = (It - 2/3*dt*Pr*Lop(:,:));
    A_trunc = A_full(2:Nc,2:Nc);
    wm3(1,m+1)=1;
    wm3(Nc+1,m+1)=0;
    B = 0;
    RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[wm3(1,m+1);wm3(Nc+1,m+1)];
    wm3(2:Nc,m+1)=A_trunc\RHS;
    % max(abs(A_full(2:Nc,:)*wm3(:,m+1)-B))
    %
    % now solve for phi
    A_full = Lop;
    A_trunc = A_full(2:Nc,2:Nc);
    phim3(1,m+1)=0;
    phim3(Nc+1,m+1)=0;
    B = -wm3(2:Nc,m+1);
    RHS = B - [A_full(2:Nc,1),A_full(2:Nc,Nc+1)]*[phim3(1,m+1);phim3(Nc+1,m+1)];
    phim3(2:Nc,m+1)=A_trunc\RHS;
    % max(abs(A_full(2:Nc,:)*phim3(:,m+1)-B))
    
    % % Peichun's computation
    % wm3_PT(1:Nc+1, m+1) = ( It - 2/3*dt*Pr*[ zeros(1,Nc+1); Lop(2:Nc, :); zeros(1,Nc+1)] )\ [1; zeros(Nc,1)];
    % phim3_PT(2:Nc, m+1) = -  Lop(2:Nc,2:Nc)\ wm3_PT(2:Nc, m+1);    %  To solve the matrix equation Ax = b, enter x=A\b
    % phim3_PT(1, m+1) = 0; phim3_PT(Nc+1, m+1) = 0;
    % comare her solution to my solution:
    % max(abs(wm3(:,m+1)-wm3_PT(:,m+1)))
    % max(abs(phim3(:,m+1)-phim3_PT(:,m+1)))
    
    % calculating the corresponding Influence Matrix M.  First row
    % is the radial derivatives of phim1 and phim3 at the inner
    % electrode.  The second row is the radial derivatives of phim1
    % and phim3 at the outer electrode.
    Mi(:,:, m+1) = [ 2*D(Nc+1,:)*phim1(:, m+1), 2*D(Nc+1, :)*phim3(:, m+1); 2*D(1,:)*phim1(:, m+1), 2*D(1,:)*phim3(:, m+1)];
    %            % Mi(:,:,m+1) = inv(Mi(:,:,m+1));   % Mi is the inverse of the influence Matrix method % Eqn. 6.72 on P.180
end

