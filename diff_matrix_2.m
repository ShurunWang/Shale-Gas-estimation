function [matrix, x] = diff_matrix_2(a,D,h1,h2,k1,k2,N,Dx)

e=ones(N+1,1);
matrix=spdiags([e -2*e e],-1:1,N+1,N+1);
x = linspace(0,Dx,N+1); % dimensionless coord
dx=Dx/N;
matrix(1,1)=h2-h1/dx;
matrix(1,2)=h1/dx;
matrix(N+1,N)=-k1/dx;
matrix(N+1,N+1)=k1/dx+k2;

for i=2:N 
    matrix(i,i-1) = -(D(i)+D(i-1))*a(i)/2/dx/dx;
    matrix(i,i) = (2*D(i)+D(i-1)+D(i+1))*a(i)/2/dx/dx;
    matrix(i,i+1) = -(D(i)+D(i+1))*a(i)/2/dx/dx;  
end