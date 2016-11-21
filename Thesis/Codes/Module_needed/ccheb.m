% function [D,DD,x] = cheb(N)
% given N, this produces x = the N+1 collocation points in [-1,1] and
% D = the (N+1)x(N+1) differentiation matrix.
% DD = the (N+1)x(N+1) matrix for the second derivative.
function [D,DD,x] = ccheb(N)

if N==0, D=0; x=1; return, end
% these are the N+1 collocation points (Peyret (3.8))
x = cos(pi*(0:N)/N)';
% these are the constants \overline{c} (Peyret (3.15))
c = [2; ones(N-1,1); 2];
% the constants are multiplied by (-1)^i in to make the formula for D
% cuter.
c = c.*(-1).^(0:N)';
% From Peyret: if for i not equal to j:
%    D(i,j) = c(i)/c(j)  (-1)^(i+j) 1/(x(i)-x(j))
% to determine D(i,i) just use the fact that the columns of D need to 
% add to zero.
%
% this matrix has X(i,j) = x(i)
X = repmat(x,1,N+1);
% this matrix has dX(i,j) = x(i)-x(j)
dX = X-X';
% the following is exactly correct for the off diagonal entries.  The
% identity matrix (eye(N+1)) is included just so that the diagonal entries
% aren't NaN.
D  = (c*(1./c)')./(dX+(eye(N+1)));     
% given the correct off diagonals, define the diagonal entries by demanding
% that the columns sum to zero.
D  = D - diag(sum(D'));                

% DD is the matrix for the second derivative.  See Peyret (3.47) for the
% formulae.  I don't have the patience to find a clever, vectorized way of
% defining DD.
c = [2; ones(N-1,1); 2];
% this is the matrix for the second derivative
i=1;
DD(i,i) = (N^4-1)/15;
for j=2:N+1
    DD(i,j)=(2/3)*((-1)^(j-1)/c(j))*((2*N^2+1)*(1-x(j))-6)/((1-x(j))^2);
end
for i=2:N
    for j=1:N+1
        if i==j
            DD(i,j) = -((N^2-1)*(1-x(i)^2)+3)/(3*(1-x(i)^2)^2);
        else
            DD(i,j) = ((-1)^(i+j)/c(j))*(x(i)^2+x(i)*x(j)-2)/((1-x(i)^2)*(x(i)-x(j))^2);
        end
    end
end
i = N+1;
for j=1:N
    DD(i,j) = (2/3)*((-1)^(j-1+N)/c(j))*((2*N^2+1)*(1+x(j))-6)/((1+x(j))^2);
end
DD(i,i) = (N^4-1)/15;


end
