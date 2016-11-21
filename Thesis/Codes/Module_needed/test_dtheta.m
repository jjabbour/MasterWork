
Nc = 24; K  = 32; q = K/2;
rs    = cos(pi*(0:Nc)/Nc).';
theta = linspace(0,2*pi,2*K+2);theta = theta(1:end-1);

% u1     = cos(rs)*sin(theta);    u1dt = cos(rs)*cos(theta);
% u2     = cos(rs)*cos(theta);    u2dt = -cos(rs)*sin(theta);        
% u3     = sin(rs)*sin(theta);    u3dt = sin(rs)*cos(theta);
% u4     = sin(rs)*cos(theta);    u4dt = -sin(rs)*sin(theta); 

u1     = cos(rs)*(exp(cos(theta)));    u1dt = -cos(rs)*((exp(cos(theta))).*sin(theta));
u2     = cos(rs)*cos(theta);    u2dt = -cos(rs)*sin(theta);        
u3     = sin(rs)*sin(theta);    u3dt = sin(rs)*cos(theta);
u4     = sin(rs)*cos(theta);    u4dt = -sin(rs)*sin(theta); 

U = mat_2_vec(u1,u2,u3,u4,Nc,K);
% Up = dtheta_electro_fft(U,K,Nc);
Up = dtheta_electro(U,K,Nc);
[u1p,u2p,u3p,u4p] = vec_2_mat(Up,Nc,K);

check1 = u1p-u1dt;  
check2 = u2p-u2dt;
check3 = u3p-u3dt;
check4 = u4p-u4dt;

disp([max(max(check1)),max(max(check2)),max(max(check3)),max(max(check4))])

